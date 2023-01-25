rm -r data/db/sequences*
python generate_files.py
./bowtie2-2.5.0-macos-arm64/bowtie2-build --threads 8 data/sequences.fasta data/db/sequences
