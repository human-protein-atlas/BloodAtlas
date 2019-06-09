<?php
#############################
## Consensus RNAseq category

log_header("\nConsensus RNAseq gene abundance and enhanced score:");

$consensus_assay = array('consensus' => 10, 'hpa' => 15, 'gtex' => 16, 'fantom' => 17, 'hpa_celline' => 19);
foreach ($consensus_assay AS $sample_type_name => $sample_type_id) {
	log_header("Gene Category calc for $sample_type_name");
	if($sample_type_name=='hpa_celline') {
		$exp_type = 'norm_exp';
	} else {
		$exp_type = 'exp';
	}
	$sql = "INSERT INTO rna_gene$temp_suffix (eg_id, sample_type_id, enhanced_score, abundance_id)
					SELECT eg_id, sample_type_id, enhanced_score,
					IF(distribution=0, '".not_detected."',
					IF(distribution=1, '".detected_in_single."',
					IF(distribution>1 AND distribution/num_tissues<0.33, '".detected_in_some."',
					IF(distribution/num_tissues>=0.33 AND distribution<num_tissues, '".detected_in_many."',
					IF(distribution=num_tissues, '".detected_in_all."', 'error!'
					))))) AS abundance_id
				FROM (
					SELECT DISTINCT eg_id, cg.sample_type_id, 
						IF(MAX($exp_type)>$lowest_nx_limit,MAX($exp_type),$lowest_nx_limit)/IF(AVG($exp_type)>$lowest_nx_limit,AVG($exp_type),$lowest_nx_limit) enhanced_score,
						COUNT(DISTINCT IF($exp_type>=$nx_detection_limit, cg.content_id, NULL)) AS distribution,
						COUNT(DISTINCT cg.content_id) AS num_tissues
					FROM consensus_gene_group$temp_suffix cg
					JOIN ensembl_genes USING (eg_id)
					".($sample_type_name!='hpa_celline' ? 'JOIN rna_category_tissue_selection rc ON rc.content_id=cg.content_id AND rc.sample_type_id=10' : '')."
					WHERE cg.sample_type_id=$sample_type_id
					GROUP BY eg_id
				) AS tmp";
	log_query("$sample_type_name RNA gene abundance and enhanced score", $sql);
	sql_query("OPTIMIZE TABLE rna_gene$temp_suffix");

	## Consensus RNAseq gene category and tissue specificity ##
	$sql = "SELECT eg_id, GROUP_CONCAT($exp_type ORDER BY $exp_type DESC), GROUP_CONCAT(cg.content_id ORDER BY $exp_type DESC), AVG($exp_type)
						FROM rna_gene$temp_suffix
						JOIN consensus_gene_group$temp_suffix cg USING (eg_id, sample_type_id)
						".($sample_type_name!='hpa_celline' ? 'JOIN rna_category_tissue_selection rc ON rc.content_id=cg.content_id AND rc.sample_type_id=10' : '')."
						WHERE cg.sample_type_id=$sample_type_id
						GROUP BY eg_id";
	$res = sql_query($sql,__LINE__);
	log_write("$sample_type_name category calc using limit: $nx_detection_limit\n");
	list($ts, $category) = get_category($res, $nx_detection_limit, $lowest_nx_limit, 5, 4);

	foreach($category as $name => $identifiers) {
		$sql = "UPDATE rna_gene$temp_suffix rg
						JOIN rna_category c ON c.category_name = '$name'
						SET rg.category_id=c.category_id 
						WHERE sample_type_id=$sample_type_id AND eg_id IN (".implode(',',$identifiers).")";
		log_query($name, $sql);
	}

	foreach($ts as $ens_id => $ts_arr) {
		$sql = "UPDATE rna_gene$temp_suffix g
						SET g.ts_score = {$ts_arr['ts_score']}
						WHERE sample_type_id=$sample_type_id AND eg_id=$ens_id";
		sql_query($sql,__LINE__);

		$insert_arr = array();
		foreach ($ts_arr['ts_tissue'] AS $content_id) {
			$insert_arr[] = "($ens_id, $sample_type_id, $content_id)";
		}
		$sql = "INSERT INTO rna_gene_".($sample_type_name=='hpa_celline' ? 'celline' : 'tissue')."_specific$temp_suffix (eg_id, sample_type_id, ".($sample_type_name=='hpa_celline' ? 'celline_id' : 'content_id').")
						VALUES ".implode(',', $insert_arr);
		sql_query($sql);
	}
}

## Functions #
function get_category($sql_res, $cutoff_detected=1, $lowest_exp_limit=0.1, $group_enriched_limit=5, $fold=5) {
	$ts = $category = array();
	while (list($ens_id, $exps, $tissues, $exp_avg) = sql_fetch($sql_res)) {
		$exp = explode(',',$exps);
		$tissue = explode(',',$tissues);
		$max = ($exp[0]>$lowest_exp_limit)?$exp[0]:$lowest_exp_limit;
		$second_max = ($exp[1]>$lowest_exp_limit)?$exp[1]:$lowest_exp_limit;
		# Not detected
		if ($max < $cutoff_detected) {
			$category['Not detected'][] = $ens_id;
			continue;
		}
		# Tissue enriched
		if (($max/$second_max) >= $fold) {
			$category['Tissue enriched'][] = $ens_id;
			$ts[$ens_id]['ts_score'] = $max/$second_max;
			$ts[$ens_id]['ts_tissue'][] = $tissue[0];
			continue;
		}
		# Group enriched
		$lim = $exp[0]/$fold;
		if ($exp[1]>=$cutoff_detected AND $exp[$group_enriched_limit]<$lim) {
			$i=1;
			while ($exp[$i]>$lim AND $i<=$group_enriched_limit) $i++;
			if ($i>1 AND $i<=$group_enriched_limit AND $exp[$i-1]>=$cutoff_detected ) {
				# check if mean(group) >= $fold
				$max_other = ($exp[$i]>=$lowest_exp_limit)?$exp[$i]:$lowest_exp_limit;
				$mean_group = array_sum(array_slice($exp,0,$i))/$i;
				if (($mean_group/$max_other) >= $fold) {
					$category['Group enriched'][] = $ens_id;
					$ts[$ens_id]['ts_score'] = $mean_group/$max_other;
					$ts[$ens_id]['ts_tissue'] = array_slice($tissue,0,$i);
					continue;
				}
			}
		}
		# Tissue enhanced
		$tissue_enhanced = array();
		foreach ($exp as $i => $fp) {
			if ($fp >= $exp_avg*$fold AND $fp >= $cutoff_detected) {
				$tissue_enhanced[] = $tissue[$i];
			}
		}
		if (!empty($tissue_enhanced)) {
			$category['Tissue enhanced'][] = $ens_id;
			$ts[$ens_id]['ts_score'] = 0;
			$ts[$ens_id]['ts_tissue'] = $tissue_enhanced;
			continue;
		}
		# Low tissue specificity (if no other case)
		$category['Low tissue specificity'][] = $ens_id;
		continue;
	}
	return array($ts, $category);
}