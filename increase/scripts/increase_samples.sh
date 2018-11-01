#!/usr/bin/env bash


increase_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

infile="$1"
analysis="$2"  # zoomin zoomout
genome="$3"  # hg19 exons tair10 saccer3
kmer=$4  # 5 3
nucleosomes_file="$5"
genome_counts="$6"
output_folder="$7"
closer_dyads_file="$8"

if [ "${analysis}" == "zoomin" ]
then
    half_window=58
else
    half_window=1000
fi


mkdir -p ${output_folder}

# Format the mutations file
python ${increase_dir}/fmt.py ${infile} ${output_folder}/mutations_unsorted.bed.gz ${genome} ${kmer}
zcat ${output_folder}/mutations_unsorted.bed.gz | sort -k1,1 -k2,2n | gzip > ${output_folder}/mutations.bed.gz

nsamples=`zcat ${output_folder}/mutations.bed.gz |cut -f4 | sort | uniq | wc -l`
nnucleosomes=`zcat ${nucleosomes_file} | wc -l`


zcat ${nucleosomes_file} | \
    awk -v half_window=${half_window} '{ OFS="\t"; }{print $1, $2-half_window, $3+half_window, $1 "_" $2 "_" $3 }' | \
    grep -v "-" | \
    intersectBed -a ${output_folder}/mutations.bed.gz -b stdin -wo -sorted | \
    awk -v half_window=${half_window} '{ OFS="\t";}{ print $9, $2-$7-half_window, $5}' | \
    gzip > ${output_folder}/mutations_dyads.tsv.gz

nmuts=`zcat ${output_folder}/mutations_dyads.tsv.gz | wc -l`
if [ ${nmuts} -lt 500 ]
then
    echo "Not enough mutations. Execution halted for ${infile}"
    exit 1
fi

zcat ${output_folder}/mutations_dyads.tsv.gz | cut -f1 | sort | uniq -c | \
    awk '{OFS ="\t"}{print $2, $1}' |  gzip > ${output_folder}/muts_per_nucl.tsv.gz


python ${increase_dir}/expected_nucleosomes.py ${genome} ${output_folder}/mutations.bed.gz \
    ${genome_counts} ${output_folder}/muts_per_nucl.tsv.gz \
    ${output_folder}/expected.npz --size ${half_window} --kmer ${kmer} --rands 1000


if [ "${analysis}" == "zoomin" ]
then
    python ${increase_dir}/increase.py --${analysis} \
        --samples ${nsamples} \
        --nucleosomes ${nnucleosomes} \
        ${output_folder}/mutations_dyads.tsv.gz \
        ${output_folder}/expected.npz \
        ${output_folder}
else
    python ${increase_dir}/increase.py --${analysis} \
        --samples ${nsamples} \
        --nucleosomes ${nnucleosomes} \
        --dyads ${closer_dyads_file} \
        ${output_folder}/mutations_dyads.tsv.gz \
        ${output_folder}/expected.npz \
        ${output_folder}
fi