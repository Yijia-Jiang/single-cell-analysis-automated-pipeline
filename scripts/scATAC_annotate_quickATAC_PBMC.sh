#!/bin/bash

spikefile=spikein_fragment_PBMC_narrowPeak.tsv # https://www.dropbox.com/scl/fi/15fg4cnk3t8ractquvm07/spikein_fragment_PBMC_narrowPeak.tsv?rlkey=jges942v31lkpt1es4pv3llyp&dl=0
bedfile=PBMC_merged_peaks.narrowPeak.filter.sorted.bed # https://www.dropbox.com/scl/fi/fbtgyv38yn24w3divzkpr/PBMC_merged_peaks.narrowPeak.filter.sorted.bed?rlkey=uze4tkc4hse6azwj70g3hsx1m&dl=0
chromsize_file=hg38.chrom.sizes.txt # https://www.dropbox.com/scl/fi/x9i1fpno5n7u9wus92l0t/hg38.chrom.sizes.txt?rlkey=3aldjof2gjp3ahm23rgui70b3&dl=0

sample=lung_SM-A62E9_1
fragFile=input/$sample/*fragments.tsv
barcodeFile=input/$sample/*fragments.tsv


do echo "start processing for sample $sample" && \
quick filter-barcodes input/$sample/*fragments.tsv --barcodes input/$sample/*barcodes.tsv > input/$sample/cell_filtered_fragments.tmp.tsv && \
cat input/$spikefile input/$sample/cell_filtered_fragments.tmp.tsv > input/$sample/concat_fake_true_fragments.tmp.tsv && \
sort -k 1,1 -k2,2n input/$sample/concat_fake_true_fragments.tmp.tsv > input/$sample/concat_fake_true_fragments.sorted.tsv && \
echo "finish sorting concated spike-in fragments" && mkdir counts/$sample'_vPBMC'/ && 
quick agg-countmatrix input/$sample/concat_fake_true_fragments.sorted.tsv -g input/$chromsize_file -p input/$bedfile --max-fragsize 999999999 -o counts/$sample'_vPBMC'/ && \
echo "Done!"
