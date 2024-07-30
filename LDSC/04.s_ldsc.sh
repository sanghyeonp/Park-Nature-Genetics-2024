#!/bin/bash

MungeName=MetS_GWAS_munge
python2 ./ldsc/ldsc.py --h2 ${MungeName}.sumstats.gz \
                    --ref-ld-chr ./baseline_v1.2/baseline. \
                    --w-ld-chr ./1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
                    --overlap-annot \
                    --frqfile-chr ./1000G_Phase3_frq/1000G.EUR.QC. \
                    --out ${MungeName}_partitioned_h2
