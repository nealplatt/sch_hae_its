#!/bin/bash

cd /master/nplatt/sch_hae_its-nigeria/data/seq_data

# Path to your accessions list
ACCESSION_LIST="accessions.list"

# Loop through each accession
while read -r ACCESSION; do
    echo "Processing $ACCESSION..."

    # Convert to FASTQ (using 4 threads, splits for paired-end)
    conda run -n sra-tools fasterq-dump "$ACCESSION" --split-files --threads 12 &

    echo "$ACCESSION done."
done < "$ACCESSION_LIST"

# parallel -j 12 gzip {} ::: *.fastq