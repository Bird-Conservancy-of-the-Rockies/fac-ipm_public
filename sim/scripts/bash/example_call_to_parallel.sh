# /bin/bash
parallel --jobs 5 'Rscript scripts/R/sim_master_effBreed.r basic effBreed_23dec2023 noZ 10 100 100 noAnn {}' ::: {1..10}

## Repeat the above call, replacing the model name argument ('effBreed_23dec2023' in the above call) with each of:
## effBreed_23dec2023
## effBreed_brPhi_23dec2023
## effBreed_brPhiNbPhi_23dec2023
## effBreed_brPhiRho_23dec2023
## effBreed_brNbPhiRho_23dec2023
## effBreed_brNbPhiRho_tight_23dec2023

