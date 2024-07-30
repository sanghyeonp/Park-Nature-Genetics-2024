#!/bin/bash
source 00_global_variables.sh

inDelimiter=tab
SelectedColumns='SNP MAF A1 A2 Effect SE Pval'
RenamedColumns='SNP MAF A1 A2 EFFECT SE P'
outDelimiter=tab
outFileDeco=.REFORMAT4ldsc.GWAS
outExtension=txt
outDIR=./reformatted/

mkdir -p ${outDIR}

for i in ${!Phenotypes[@]}; do
    phenotype=${Phenotypes[$i]}
    QCsumstat=${Cleaned_sumstat[$i]}

    echo "Reformatting ${phenotype}"

    python3 exe.reformat_sumstat.py --sumstat ${QCsumstat}\
                            --in_delimiter ${inDelimiter}\
                            --col_selection ${SelectedColumns}\
                            --col_rename ${RenamedColumns}\
                            --out_delimiter ${outDelimiter}\
                            --suffix ${outFileDeco}\
                            --out_file_extension ${outExtension}\
                            --outdir ${outDIR}\
                            --verbose
done


## For no UKB
inDelimiter=tab
SelectedColumns='SNP MAF A1 A2 Effect SE Pval'
RenamedColumns='SNP MAF A1 A2 EFFECT SE P'
outDelimiter=tab
outFileDeco=.REFORMAT4ldsc.noUKB.GWAS
outExtension=txt
outDIR=./reformatted/


mkdir -p ${outDIR}

for i in ${!Phenotypes[@]}; do
    phenotype=${Phenotypes[$i]}
    QCsumstat=${Cleaned_sumstat_noUKB[$i]}

    echo "Reformatting ${phenotype}"

    python3 exe.reformat_sumstat.py --sumstat ${QCsumstat}\
                            --in_delimiter ${inDelimiter}\
                            --col_selection ${SelectedColumns}\
                            --col_rename ${RenamedColumns}\
                            --out_delimiter ${outDelimiter}\
                            --suffix ${outFileDeco}\
                            --out_file_extension ${outExtension}\
                            --outdir ${outDIR}\
                            --verbose
done