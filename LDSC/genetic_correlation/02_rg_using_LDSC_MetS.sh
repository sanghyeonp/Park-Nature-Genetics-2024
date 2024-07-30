#!/bin/bash

metadata=./metadata_for_LDSC_rg.csv
outDIR=./rg_ldsc

mkdir -p $outDIR


while IFS=, read -r Category Trait MetSMunge TraitMunge
do
    echo $Trait
    if [ "$Category" == "Category" ]; then
        continue
    else
        echo "LDSC rg calculation: ${Category} ${Trait}"
        Trait_reformat=${Trait//[ ]/_}
        Category_reformat=${Category//[ ]/_}

        ## MetS genetic correlation
        echo "LDSC_rg_MetS_Category_${Category_reformat}_Trait_${Trait_reformat}"
        python2 ${src_LDSC} --rg "${MetSMunge}","${TraitMunge}"\
                        --ref-ld-chr ../eur_w_ld_chr/\
                        --w-ld-chr ../eur_w_ld_chr/\
                        --out "${outDIR}/LDSC_rg_MetS_Category_${Category_reformat}_Trait_${Trait_reformat}"
    fi
done < ${metadata}

