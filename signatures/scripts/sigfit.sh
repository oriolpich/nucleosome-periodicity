#!/usr/bin/env bash

# Run sigfit over a set of mutation files and
# assign to each mutation the most probable signature
# The process is:
# - preprocess the mutations data to generate the format required by *sigfit*
# - use *sigfit* to compute the weights of each sample for each signature
# - compute the most probable signature contributing to each mutation


this_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

source activate env_nucperiod

mutations_folder="$1"
project_folder="$2"

mkdir -p ${project_folder}

# This part can be parallelized
for mut_file in ${mutations_folder}/*.tsv.gz
do

    name=$(basename ${mut_file})
    name=${name/.tsv.gz/}
    output_folder=${project_folder}/${name}
    mkdir -p ${output_folder}

    sigfit_infile=${output_folder}/sigfit.tsv.gz
    sigfit_outfile=${output_folder}/exposures.txt

    # Prepare the file
    python ${this_dir}/fmt.py ${mut_file} ${sigfit_infile}
    ret=$?
    if [ $ret -eq 3 ]
    then
        echo "File ${mut_file} does not have enough mutations in any sample"
        continue
    fi

    # Run analysis
    source deactivate
    source activate env_sigfit
    Rscript ${this_dir}/extractSig.R ${sigfit_infile} ${sigfit_outfile}
    source deactivate
    source activate env_nucperiod

    # Assing mutations
    python ${this_dir}/assign.py ${sigfit_infile} ${sigfit_outfile} signatures_probabilities.txt ${output_folder}

done
