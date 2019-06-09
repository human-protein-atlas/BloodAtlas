## cutoff
under_limit <- function(new.expression, expression, method, scale_by = 1, limit_method = "GTEx") {
	tibble(new.expression, expression, method) %>%
		filter(method == limit_method) %$%
		sort(new.expression)[{a <- sort(c(expression, 1)); which(a == 1)[1]}] * scale_by
}

## KNN: missing data imputation
impute_expression <- function(expression, tissue.method, ensg_id) {
  # Impute missing values
	write("  -- Impute missing values", stdout())
  tibble(tissue.method, expression, ensg_id) %>%
    left_join(spread(., key = tissue.method, value = expression) %>% 
                column_to_rownames("ensg_id") %>% 
                as.matrix() %>% 
                impute::impute.knn() %>% 
                {data <- .$data; ensg_id. <- rownames(data); as.tibble(data) %>% mutate(ensg_id = ensg_id.)} %>%
                gather(key = "tissue.method", value = "imputed.expression", -ensg_id),
              by = c("ensg_id", "tissue.method")) %$%
    return(imputed.expression)
}

# TMM: trimmed mean of M-values normalization within methods
tmm_method_normalization <- function(expression, method, tissue.method, ensg_id){
  methods <- unique(method)
  tibble(expression, tissue.method, ensg_id) %>%
    left_join({for(i in 1:length(methods)){
      tb <-
        filter(., method == methods[i]) %>%
        spread(key = tissue.method, value = expression) %>%
        column_to_rownames("ensg_id") %>%
        as.matrix() %>%
        NOISeq::tmm() %>%
        {names <- rownames(.); as.tibble(.) %>% mutate(ensg_id = names)} 
      
      if(i == 1) {
        MM <- tb
      } else MM <- full_join(MM, tb, by = "ensg_id")
    }
      MM %>%
        gather(key = "tissue.method", value = "tmm.expression", -ensg_id)}, 
    by = c("ensg_id", "tissue.method")) %>% 
    mutate(tmm.expression = ifelse(is.na(tmm.expression), 
                                   expression, 
                                   tmm.expression),
           tmm.expression = ifelse(tmm.expression<0, 0, 
                                   tmm.expression)) %$%
    return(tmm.expression)
}

#SD: Method and gene standard deviation on dataset scaled values
pareto_scale_method_gene <- function(expression, method, ensg_id) {
  tibble(method, expression, ensg_id) %>%
    left_join(group_by(., method, ensg_id) %>%
                summarise(method.ensg_id.sd = sd(expression, na.rm = T)), 
              by = c("method", "ensg_id")) %>%
    
    # PAR: Method and gene pareto scaling on dataset scaled values
    mutate(norm.expression = expression/sqrt(method.ensg_id.sd),
           norm.expression = ifelse(!is.finite(norm.expression), 0, norm.expression)) %$%
    return(norm.expression)
}

# LIMMA: Batch correction of method, method and gene pareto scaled values without blood
limma_method_correction <- function(expression, method, tissue.method, ensg_id, filtered.methods){
  tibble(expression, tissue.method, ensg_id) %>%
    left_join(filter(., !method %in% filtered.methods) %>%
                spread(key = tissue.method, value = expression) %>%
                column_to_rownames("ensg_id") %>%
                as.matrix() %>%
                {.[apply(., 1, function(x) !any(is.na(x))),]} %>%
                {log(.+1, 10)} %>%
                limma::removeBatchEffect(batch = as.factor(str_extract(colnames(.), "(?<=\\.).*"))) %>%
                {names <- rownames(.); as.tibble(.) %>% mutate(ensg_id = names)} %>%
                gather(key = "tissue.method", value = "method.corrected.expression", -ensg_id) %>%
                mutate(method.corrected.expression = 10^method.corrected.expression - 1,
                       method.corrected.expression = ifelse(is.na(method.corrected.expression),
                                                            norm.dsscaled.expression, 
                                                            method.corrected.expression),
                       method.corrected.expression = ifelse(method.corrected.expression<0, 0,
                                                            method.corrected.expression)), 
              by = c("ensg_id", "tissue.method")) %>% 
    mutate(method.corrected.expression = ifelse(method %in% filtered.methods, 
                                                expression, 
                                                method.corrected.expression)) %$%
    return(method.corrected.expression)
}