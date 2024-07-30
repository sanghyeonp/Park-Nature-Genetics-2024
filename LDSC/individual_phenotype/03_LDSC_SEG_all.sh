#!/bin/bash
source 00_global_variables.sh
saveto_parent=./LDSC_SEG_analysis
mkdir -p $saveto_parent


cts_Cahoy=../Cahoy.ldcts
cts_Multi_tissue_gene_expr=../Multi_tissue_gene_expr.ldcts
cts_Multi_tissue_chromatin=../Multi_tissue_chromatin.ldcts
cts_GTEx_brain=../GTEx_brain.ldcts
cts_ImmGen=../ImmGen.ldcts
cts_Corces_ATAC=../Corces_ATAC.ldcts
LDSCSEG_weights=../weights_hm3_no_hla/weights.

for idx1 in ${!Phenotypes[@]}; do
    phenotype=${Phenotypes[$idx1]}

    saveto="${saveto_parent}/${phenotype}"
    mkdir -p $saveto

    Munged_sumstat=${munged_sumstat[$idx1]}

    cts_list=('Cahoy' 'Multi_tissue_gene_expr' 'Multi_tissue_chromatin' 'GTEx_brain')
    cts_ld_list=(${cts_Cahoy} ${cts_Multi_tissue_gene_expr} ${cts_Multi_tissue_chromatin}
            ${cts_GTEx_brain})

    for idx2 in "${!cts_list[@]}"; do
        python2 ../ldsc/ldsc.py --h2-cts ${Munged_sumstat} \
                            --ref-ld-chr ../1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD. \
                            --out ${saveto}/${phenotype}_${cts_list[$idx2]} \
                            --ref-ld-chr-cts ${cts_ld_list[$idx2]} \
                            --w-ld-chr ../1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
    done
done
