manhattan.data = dbGetQuery(
	db,
	paste(
		sep =' ',
		'SELECT T1.chromosome AS `chromosome`, T1.position as `position`, T2.rsid AS `rsid`, T2.alleleA AS `alleleA`, T2.alleleB AS `alleleB`, `add:meta_beta`, `add:meta_se`, `add:meta_pvalue`, T1.mean_bf AS `mean_bf`,',
		'`Gambia:type` AS "Gambia:impute_type", `Gambia:info` AS "Gambia:impute_info",',
		'( `Gambia:controls_AB` + 2 * `Gambia:controls_AA` )  AS "Gambia:controls_alleleA_count",',
		'( `Gambia:controls_AB` + 2 * `Gambia:controls_BB` )  AS "Gambia:controls_alleleB_count",',
		'( `Gambia:controls_AB` + 2 * `Gambia:controls_BB` ) / ( 2 * ( `Gambia:controls_AA` + `Gambia:controls_AB` + `Gambia:controls_BB` )) AS "Gambia:controls_alleleB_frequency",',
		'`Malawi:type` AS "Malawi:impute_type", `Malawi:info` AS "Malawi:impute_info",',
		'( `Malawi:controls_AB` + 2 * `Malawi:controls_BB` ) / ( 2 * ( `Malawi:controls_AA` + `Malawi:controls_AB` + `Malawi:controls_BB` )) AS "Malawi:controls_alleleB_frequency",',
		'`Kenya:type` AS "Kenya:impute_type", `Kenya:info` AS "Kenya:impute_info",',
		'( `Kenya:controls_AB` + 2 * `Kenya:controls_AA` )  AS "Kenya:controls_alleleA_count",',
		'( `Kenya:controls_AB` + 2 * `Kenya:controls_BB` )  AS "Kenya:controls_alleleB_count",',
		'( `Kenya:controls_AB` + 2 * `Kenya:controls_BB` ) / ( 2 * ( `Kenya:controls_AA` + `Kenya:controls_AB` + `Kenya:controls_BB` )) AS 
		"Kenya:controls_alleleB_frequency"',
		'FROM MeanBayesFactorBaked T1 INNER JOIN CrossModelCountFilteredView T2 ON T1.variant_id == T2.variant_id WHERE T1.prior_id == 1'
#		,' LIMIT 1000'
	)
)
