#!/bin/bash

NAME=MetS
LDSCOUT=../../ldscout/MetS_ldscout.RData
outDir=.
#######################################

Rscript exc.parallel_analysis.R --name "${NAME}" \
                    --rdata "${LDSCOUT}" \
                    --outDir "${outDir}"


NAME=MetS_noUKB
LDSCOUT=../../ldscout/MetS_noUKB_ldscout.RData
outDir=.
#######################################

Rscript exc.parallel_analysis.R --name "${NAME}" \
                    --rdata "${LDSCOUT}" \
                    --outDir "${outDir}"