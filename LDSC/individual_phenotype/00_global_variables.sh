Phenotypes=('BMI'
			'WC'
			'T2D'
			'FG'
			'HTN'
			'TG'
			'HDL')

Cleaned_sumstat=('CLEANED.BMI_GIANT2018.txt'
				'CLEANED.WC_UKB.txt'
				'CLEANED.T2D_meta_FinngenR7_Mahajan2022_MVP.txt'
				'CLEANED.FG_MAGIC.txt'
				'CLEANED.HTN_meta_FinngenR7_UKB.txt'
				'CLEANED.TG_GLGC_UKB.txt'
				'CLEANED.HDL_GLGC_UKB.txt'
				)

SampleSize=(806834
			385932
			597437
			151188
			508612
			1253277
			1244580
			)

reformated_sumstat_for_LDSC=('./reformatted/CLEANED.BMI_GIANT2018.REFORMAT4ldsc.GWAS.txt'
							'./reformatted/CLEANED.WC_UKB.REFORMAT4ldsc.GWAS.txt'
							'./reformatted/CLEANED.T2D_meta_FinngenR7_Mahajan2022_MVP.REFORMAT4ldsc.GWAS.txt'
							'./reformatted/CLEANED.FG_MAGIC.REFORMAT4ldsc.GWAS.txt'
							'./reformatted/CLEANED.HTN_meta_FinngenR7_UKB.REFORMAT4ldsc.GWAS.txt'
							'./reformatted/CLEANED.TG_GLGC_UKB.REFORMAT4ldsc.GWAS.txt'
							'./reformatted/CLEANED.HDL_GLGC_UKB.REFORMAT4ldsc.GWAS.txt'
							)

munged_sumstat=('./munged/BMI_GWAS_munge.sumstats.gz'
				'./munged/WC_GWAS_munge.sumstats.gz'
				'./munged/T2D_GWAS_munge.sumstats.gz'
				'./munged/FG_GWAS_munge.sumstats.gz'
				'./munged/HTN_GWAS_munge.sumstats.gz'
				'./munged/TG_GWAS_munge.sumstats.gz'
				'./munged/HDL_GWAS_munge.sumstats.gz'
				)

#############################################################################################################################################################################################

Cleaned_sumstat_noUKB=('CLEANED.BMI_GIANT2015.txt'
						'CLEANED.WC_GIANT2015.txt'
						'CLEANED.T2D_meta_FinngenR7_MVP.txt'
						'CLEANED.FG_MAGIC.txt'
						'CLEANED.HTN_Finngen_r7.txt'
						'CLEANED.TG_GLGC_noUKB.txt'
						'CLEANED.HDL_GLGC_noUKB.txt'
						)

SampleSize_noUKB=(322154
				232101
				345698
				151188
				247289
				864240
				888227
				)

reformated_sumstat_for_LDSC_noUKB=('./reformatted/CLEANED.BMI_GIANT2015.REFORMAT4ldsc.noUKB.GWAS.txt'
									'./reformatted/CLEANED.WC_GIANT2015.REFORMAT4ldsc.noUKB.GWAS.txt'
									'./reformatted/CLEANED.T2D_meta_FinngenR7_MVP.REFORMAT4ldsc.noUKB.GWAS.txt'
									'./reformatted/CLEANED.FG_MAGIC.REFORMAT4ldsc.noUKB.GWAS.txt'
									'./reformatted/CLEANED.HTN_Finngen_r7.REFORMAT4ldsc.noUKB.GWAS.txt'
									'./reformatted/CLEANED.TG_GLGC_noUKB.REFORMAT4ldsc.noUKB.GWAS.txt'
									'./reformatted/CLEANED.HDL_GLGC_noUKB.REFORMAT4ldsc.noUKB.GWAS.txt'
									)

munged_sumstat_noUKB=('./munged/BMI_GWAS_noUKB_munge.sumstats.gz'
				'./munged/WC_GWAS_noUKB_munge.sumstats.gz'
				'./munged/T2D_GWAS_noUKB_munge.sumstats.gz'
				'./munged/FG_GWAS_noUKB_munge.sumstats.gz'
				'./munged/HTN_GWAS_noUKB_munge.sumstats.gz'
				'./munged/TG_GWAS_noUKB_munge.sumstats.gz'
				'./munged/HDL_GWAS_noUKB_munge.sumstats.gz'
				)
