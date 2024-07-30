#!/bin/bash
source 00_global_variables.sh


OUTdir_for_munge=./munged
mkdir -p $OUTdir_for_munge

OUTdir_for_heritability=./heritability
mkdir -p $OUTdir_for_heritability

for i in ${!Phenotypes[@]}; do
    phenotype=${Phenotypes[$i]}

    SUMSTAT=${reformated_sumstat_for_LDSC[$i]}
    EffN=${SampleSize[$i]}
    MungeName="${phenotype}_GWAS_munge"
    LDSCName="${phenotype}_GWAS_ldsc"

    if [[ ! -f "${OUTdir_for_munge}/${MungeName}_munge.sumstats.gz" ]]; then
        python2 ../ldsc/munge_sumstats.py --sumstats ${SUMSTAT} \
                                --N ${EffN} \
                                --out ${OUTdir_for_munge}/${MungeName} \
                                --merge-alleles ../w_hm3.snplist \
                                --chunksize 500000
    fi

    python2 ../ldsc/ldsc.py --h2 ${OUTdir_for_munge}/${MungeName}.sumstats.gz \
                    --ref-ld-chr ../1000G_Phase3_ldscores/LDscore. \
                    --w-ld-chr ../1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
                    --out "${OUTdir_for_heritability}/${LDSCName}"
done



OUTdir_for_munge=./munged
mkdir -p $OUTdir_for_munge

OUTdir_for_heritability=./heritability
mkdir -p $OUTdir_for_heritability


for i in ${!Phenotypes[@]}; do
    phenotype=${Phenotypes[$i]}

    SUMSTAT=${reformated_sumstat_for_LDSC_noUKB[$i]}
    EffN=${SampleSize_noUKB[$i]}
    MungeName="${phenotype}_GWAS_noUKB_munge"
    LDSCName="${phenotype}_GWAS_noUKB_ldsc"

    if [[ ! -f "${OUTdir_for_munge}/${MungeName}_munge.sumstats.gz" ]]; then
        python2 ${src_Munge} --sumstats ${SUMSTAT} \
                                --N ${EffN} \
                                --out ${OUTdir_for_munge}/${MungeName} \
                                --merge-alleles ../w_hm3.snplist \
                                --chunksize 500000
    fi

    python2 ${src_LDSC} --h2 ${OUTdir_for_munge}/${MungeName}.sumstats.gz \
                        --ref-ld-chr ../1000G_Phase3_ldscores/LDscore. \
                        --w-ld-chr ../1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
                        --out "${OUTdir_for_heritability}/${LDSCName}"
done
