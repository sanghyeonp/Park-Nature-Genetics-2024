#!/bin/bash
source 00_global_variables.sh


cts_Cahoy=./Cahoy.ldcts
cts_Multi_tissue_gene_expr=./Multi_tissue_gene_expr.ldcts
cts_Multi_tissue_chromatin=./Multi_tissue_chromatin.ldcts
cts_GTEx_brain=./GTEx_brain.ldcts
cts_ImmGen=./ImmGen.ldcts
cts_Corces_ATAC=./Corces_ATAC.ldcts
LDSCSEG_weights=./weights_hm3_no_hla/weights.


cts_list=('Cahoy' 'Multi_tissue_gene_expr' 'Multi_tissue_chromatin' 'GTEx_brain')
cts_ld_list=(${cts_Cahoy} ${cts_Multi_tissue_gene_expr} ${cts_Multi_tissue_chromatin}
        ${cts_GTEx_brain})

MungeName=MetS_GWAS_munge

for idx2 in "${!cts_list[@]}"; do
    python2 ./ldsc/ldsc.py --h2-cts ${MungeName}.sumstats.gz \
                        --ref-ld-chr ./1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD. \
                        --out ${MungeName}_${cts_list[$idx2]} \
                        --ref-ld-chr-cts ${cts_ld_list[$idx2]} \
                        --w-ld-chr ./1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
done
