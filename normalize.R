#---
#title: "HPA Atlas normalization"
#author: "Max Karlsson"
#created date: "2018 August 28"
#modified by: "Wen Zhong"
#modified date: "2018 Nov 28"
#modified by:	"Per Oksvold"
#modified date: "2018 Dec 18"
#
# cmd: Rscript normalize.R
# dependencies: tsv files in /tmp/ folder
#	output: 
#			/tmp/norm_expression.tsv (includes both nx with and without limma)
#			/tmp/norm_pig_expression.tsv
#			/tmp/norm_mouse_expression.tsv
#			/tmp/norm_celline_expression.tsv
#---

write(paste("\nNormalize TPM R-script started (", Sys.time(), ")\n"), stdout())

suppressMessages(library('impute', quietly = TRUE))
suppressMessages(library('tidyverse', quietly = TRUE))
suppressMessages(library('limma', quietly = TRUE))
suppressMessages(library('NOISeq', quietly = TRUE))
suppressMessages(library('magrittr', quietly = TRUE))

source('function.R')


result_folder <- '/tmp/'
geneinfo_path <- '/tmp/gene_info.tsv'
hpa_path <- '/tmp/rna_hpa.tsv'
gtex_path <- '/tmp/rna_gtex.tsv'
fantom_path <- '/tmp/rna_fantom.tsv'
blood_path <- '/tmp/rna_blood.tsv'
celline_path <- '/tmp/rna_celline.tsv'

####################
## Read Data
write("  -- Read data files", stdout())

write("  -- Read gene data", stdout())
ensg.info <- 
  geneinfo_path %>%
  readr::read_delim(delim = "\t", 
		col_types = cols(
			ensg_id = col_character()
		))

write("  -- Read hpa tissue data", stdout())
hpa.atlas <-
  hpa_path %>%
  readr::read_delim(delim = "\t", 
		col_types = cols(
			ensg_id = col_character(),
			tissue = col_character(),
			expression = col_double(),
			sample_type_id = col_integer()
		)) %>%
  dplyr::rename(tissue = 2,
                expression = 3)

write("  -- Read fantom data", stdout())
fantom.atlas <-
  fantom_path %>%
  readr::read_delim(delim = "\t", 
		col_types = cols(
			ensg_id = col_character(),
			tissue = col_character(),
			expression = col_double(),
			sample_type_id = col_integer()
		)) %>%
  dplyr::rename(tissue = 2,
                expression = 3) %>%
  right_join(crossing(ensg_id = ensg.info$ensg_id, tissue = unique(.$tissue)), 
             by = c("ensg_id", "tissue"))

write("  -- Read gtex data", stdout())
gtex.atlas <-
  gtex_path %>%
  readr::read_delim(delim = "\t",
			col_types = cols(
			ensg_id = col_character(),
			tissue = col_character(),
			expression = col_double(),
			sample_type_id = col_integer()
		)) %>%
  dplyr::rename(tissue = 2,
                expression = 3) %>%
  right_join(crossing(ensg_id = ensg.info$ensg_id, tissue = unique(.$tissue)), 
             by = c("ensg_id", "tissue"))

write("  -- Read blood data", stdout())
blood.atlas <-
  blood_path %>%
  readr::read_delim(delim = "\t", 
			col_types = cols(
			ensg_id = col_character(),
			tissue = col_character(),
			expression = col_double(),
			sample_type_id = col_integer()
		))
				 
write("  -- Read hpa celline data", stdout())
celline.raw <-
  celline_path %>%
  readr::read_delim(delim = "\t", 
		col_types = cols(
			ensg_id = col_character(),
			tissue = col_character(),
			expression = col_double(),
			sample_type_id = col_integer()
		)) %>%
  mutate(method = factor(c(rep("CELLINE")),
										levels = c("CELLINE"))) %>%
  mutate(expression = round(expression, 4), 
         tissue.method = paste(tissue, method, sep = "."))

## combine atlas datasets
write("-- Combine consensus tissue atlas datasets", stdout())
all.atlas.raw <- 
  rbind(dplyr::select(hpa.atlas, ensg_id, tissue, expression, sample_type_id), 
        dplyr::select(gtex.atlas, ensg_id, tissue, expression, sample_type_id), 
        dplyr::select(fantom.atlas, ensg_id, tissue, expression, sample_type_id),
        dplyr::select(blood.atlas, ensg_id, tissue, expression, sample_type_id)) %>%
  mutate(method = factor(c(rep("HPA", nrow(hpa.atlas)), 
                           rep("GTEx", nrow(gtex.atlas)), 
                           rep("FANTOM", nrow(fantom.atlas)),
                           rep("Blood", nrow(blood.atlas))),
                         levels = c("Blood", "HPA", "GTEx", "FANTOM"))) %>%
  mutate(expression = round(expression, 4), 
         tissue.method = paste(tissue, method, sep = ".")) 

# Dev - show result for debugging
#write("  ---- Dev: Write raw data tissue file", stdout())
#readr::write_delim(all.atlas.raw, path = paste0(result_folder, "all.atlas.raw.tsv"), delim = "\t")


################
## Normalization

## Tissue consensus normalization
write("-- Tissue Consensus Normalization", stdout())
all.atlas <- 
  all.atlas.raw %>%
  mutate(
    # Impute missing values
    imputed = case_when(is.na(expression) ~ TRUE,
                        TRUE ~ FALSE),
		imputed.zero.expression = ifelse(imputed, 0, expression),
    # TMM scaling of data with imputation (set to 0)
    dstmm.zero.expression = tmm_method_normalization(imputed.zero.expression, method, tissue.method, ensg_id),
		# Gene pareto. "Imputed" values are set to NA to not count.
		gene_dstmm.zero.impute.expression = pareto_scale_method_gene(ifelse(imputed, NA, dstmm.zero.expression), 
																																 method, ensg_id),
		# Limma
		limma_gene_dstmm.zero.impute.expression = limma_method_correction(gene_dstmm.zero.impute.expression, method,
																																			tissue.method, ensg_id,
																																			filtered.methods = "Blood"))  %>%
    # Scale so that under limit is 1
    mutate_at(.funs = funs(. / (under_limit(., expression, method, scale_by = 1))), 
              .vars = grep(".expression$", colnames(.), value = T)) 

# Dev - show result for debugging
#write("  ---- Dev: Write normalized tissue file", stdout())
#readr::write_delim(all.atlas, path = paste0(result_folder, "all.atlas.tsv"), delim = "\t")

# Only select relevant columns for different datasets. Get both NX values: with and without limma
write("  -- Get clean dataset", stdout())
all.atlas.clean <- select(all.atlas, ensg_id, tissue, sample_type_id, gene_dstmm.zero.impute.expression, limma_gene_dstmm.zero.impute.expression)

# Write result to file
write("  -- Write result to file", stdout())
readr::write_delim(all.atlas.clean, path = paste0(result_folder, "norm_expression.tsv"), delim = "\t")


## Methods to normalize without limma
run_methods <- c('celline')
temp_data <- data.frame();
for (current_method in run_methods) {
    ## normalization without limma
		write(paste0("-- ", current_method," Normalization"), stdout())
		## make sure temp_data does not contain old values
		rm(temp_data)
		temp_data <- 
			get(paste0(current_method, '.raw')) %>%
			mutate(
				# TMM scaling of data
				dstmm.expression = tmm_method_normalization(expression, method, tissue.method, ensg_id),
				# Gene pareto
				gene_dstmm.expression = pareto_scale_method_gene(dstmm.expression, method, ensg_id))  %>%
				# Scale so that under limit is 1, use current method (only one in each dataframe)
				mutate_at(.funs = funs(. / (under_limit(., expression, method, scale_by = 1, limit_method = toupper(current_method)))), 
									.vars = grep(".expression$", colnames(.), value = T)) 

		assign(paste0(current_method, '.norm'), temp_data)

		# Dev - show result for debugging
		#write("  ---- Dev: Write normalized file", stdout())
		#readr::write_delim(get(paste0(current_method, '.norm')), path = paste0(result_folder, current_method, ".norm.tsv"), delim = "\t")
		
		# Only select relevant columns for different datasets
		write("  -- Get clean dataset", stdout())
		assign(paste0(current_method, '.norm.clean'), select(get(paste0(current_method, '.norm')), ensg_id, tissue, sample_type_id, gene_dstmm.expression))

		# Write result to file
		write("  -- Write result to file", stdout())
		readr::write_delim(get(paste0(current_method, '.norm.clean')), path = paste0(result_folder, "norm_", current_method, "_expression.tsv"), delim = "\t")
}

write(paste("\nFinished (", Sys.time(), ")\n"), stdout())
