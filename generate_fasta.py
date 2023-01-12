import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import progressbar


n = 1000000
csv.field_size_limit(100000000)

try:
    os.remove("sequences.fasta")
except FileNotFoundError:
    pass

with progressbar.ProgressBar(max_value=n) as bar:
    with open("sequences.fasta", "w") as fasta_file:
        with open("sequences.csv") as csv_file:
            reader = csv.DictReader(csv_file)
            count = 0
            for row in reader:
                record = SeqRecord(
                    Seq(row["sequence"]),
                    id=row["id"],
                    description=row["target_gene"]
                )
                SeqIO.write([record], fasta_file, "fasta")
                count = count + 1
                bar.update(count)
                if count >= n:
                    break
