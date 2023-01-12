rm -r db/sequences*
./bowtie2-2.5.0-macos-arm64/bowtie2-build --threads 8 sequences.fasta db/sequences
#./bowtie2-2.5.0-macos-arm64/bowtie2-inspect --summary db/sequences

