#!/bin/bash

NAME=MetS_noUKB
LDSCOUT=../../2_ldsc/ldscout/MetS_noUKB_ldscout.RData
outDir=.
#######################################

Rscript exc.compute_KMO.R --name "${NAME}" \
                    --rdata "${LDSCOUT}" \
                    --outDir "${outDir}"

NAME=MetS
LDSCOUT=../../2_ldsc/ldscout/MetS_ldscout.RData
outDir=.
#######################################

Rscript exc.compute_KMO.R --name "${NAME}" \
                    --rdata "${LDSCOUT}" \
                    --outDir "${outDir}"