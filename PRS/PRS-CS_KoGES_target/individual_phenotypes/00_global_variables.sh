#!/bin/bash

Phenotypes=('BMI'
			'WC'
			'T2D'
			'FG'
			'HTN'
			'TG'
			'HDL')

Cleaned_sumstat_noUKB=('CLEANED.BMI_GIANT2015.txt'
						'CLEANED.WC_GIANT2015.txt'
						'CLEANED.T2D_meta_FinngenR7_MVP.txt'
						'CLEANED.FG_MAGIC.txt'
						'CLEANED.HTN_Finngen_r7.txt'
						'CLEANED.TG_GLGC_noUKB.txt'
						'CLEANED.HDL_GLGC_noUKB.txt'
						)

Sumstat_for_PRS_noUKB=('CLEANED.BMI_GIANT2015.REFORMAT4PRS.GWAS.txt'
						'CLEANED.WC_GIANT2015.REFORMAT4PRS.GWAS.txt'
						'CLEANED.T2D_meta_FinngenR7_MVP.REFORMAT4PRS.GWAS.txt'
						'CLEANED.FG_MAGIC.REFORMAT4PRS.GWAS.txt'
						'CLEANED.HTN_Finngen_r7.REFORMAT4PRS.GWAS.txt'
						'CLEANED.TG_GLGC_noUKB.REFORMAT4PRS.GWAS.txt'
						'CLEANED.HDL_GLGC_noUKB.REFORMAT4PRS.GWAS.txt'
						)

SampleSize_noUKB=(322154
				232101
				345698
				151188
				247289
				864240
				888227
				)
