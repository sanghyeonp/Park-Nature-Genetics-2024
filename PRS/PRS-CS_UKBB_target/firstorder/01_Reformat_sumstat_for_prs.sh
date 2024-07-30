#!/bin/bash

########################################################################################################
QCsumstat=MetS_noUKB_F1_GWAS.csv
inDelimiter=comma
SelectedColumns='SNP A1 A2 est SE Pval_Estimate'
RenamedColumns='MarkerName A1 A2 Beta SE P'
outDelimiter=tab
Suffix=.REFORMAT4PRS.GWAS
outExtension=txt
outDIR=.
########################################################################################################

python3 ../exe.reformat_sumstat.py --sumstat ${QCsumstat}\
                        --in_delimiter ${inDelimiter}\
                        --col_selection ${SelectedColumns}\
                        --col_rename ${RenamedColumns}\
                        --out_delimiter ${outDelimiter}\
                        --suffix ${Suffix}\
                        --out_file_extension ${outExtension}\
                        --outdir ${outDIR}\
                        --verbose

########################################################################################################
QCsumstat=MetS_noUKB_F2_GWAS.csv
inDelimiter=comma
SelectedColumns='SNP A1 A2 est SE Pval_Estimate'
RenamedColumns='MarkerName A1 A2 Beta SE P'
outDelimiter=tab
Suffix=.REFORMAT4PRS.GWAS
outExtension=txt
outDIR=.
########################################################################################################

python3 ../exe.reformat_sumstat.py --sumstat ${QCsumstat}\
                        --in_delimiter ${inDelimiter}\
                        --col_selection ${SelectedColumns}\
                        --col_rename ${RenamedColumns}\
                        --out_delimiter ${outDelimiter}\
                        --suffix ${Suffix}\
                        --out_file_extension ${outExtension}\
                        --outdir ${outDIR}\
                        --verbose

########################################################################################################
QCsumstat=MetS_noUKB_F3_GWAS.csv
inDelimiter=comma
SelectedColumns='SNP A1 A2 est SE Pval_Estimate'
RenamedColumns='MarkerName A1 A2 Beta SE P'
outDelimiter=tab
Suffix=.REFORMAT4PRS.GWAS
outExtension=txt
outDIR=.
########################################################################################################

python3 ../exe.reformat_sumstat.py --sumstat ${QCsumstat}\
                        --in_delimiter ${inDelimiter}\
                        --col_selection ${SelectedColumns}\
                        --col_rename ${RenamedColumns}\
                        --out_delimiter ${outDelimiter}\
                        --suffix ${Suffix}\
                        --out_file_extension ${outExtension}\
                        --outdir ${outDIR}\
                        --verbose
