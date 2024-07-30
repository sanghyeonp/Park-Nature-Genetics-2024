#!/bin/bash

source 00_global_variables.sh


OUT=.

for idx in "${!Phenotypes[@]}"
do
    sumstatPath="${Cleaned_sumstat_noUKB[$idx]}"
    inDelimiter=tab
    SelectedColumns='SNP A1 A2 Effect SE Pval'
    RenamedColumns='MarkerName A1 A2 Beta SE P'
    outDelimiter=tab
    Suffix=.REFORMAT4PRS.GWAS
    outExtension=txt

    python3 ../exe.reformat_sumstat.py --sumstat ${sumstatPath}\
                            --in_delimiter ${inDelimiter}\
                            --col_selection ${SelectedColumns}\
                            --col_rename ${RenamedColumns}\
                            --out_delimiter ${outDelimiter}\
                            --suffix ${Suffix}\
                            --out_file_extension ${outExtension}\
                            --outdir ${OUT}\
                            --verbose
done
