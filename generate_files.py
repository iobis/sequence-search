import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import progressbar
import sqlite3


# fasta file

csv.field_size_limit(100000000)

try:
    os.remove("data/sequences.fasta")
except FileNotFoundError:
    pass

with progressbar.ProgressBar(max_value=progressbar.UnknownLength) as bar:
    with open("data/sequences.fasta", "w") as fasta_file:
        with open("data/sequences.csv") as csv_file:
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

# sqlite (occurrence IDs)

try:
    os.remove("data/occurrence.sqlite")
except FileNotFoundError:
    pass

con = sqlite3.connect("data/occurrence.sqlite")
cur = con.cursor()
cur.execute("create table occurrence (id, seq)")

with progressbar.ProgressBar(max_value=progressbar.UnknownLength) as bar:
    with open("data/sequences.csv") as csv_file:
        reader = csv.DictReader(csv_file)
        count = 0
        for row in reader:
            seq = row["id"]
            ids = row["occurrence_ids"].split("|")
            cur.executemany("insert into occurrence values (?, ?)", list(zip(ids, [seq])))
            count = count + 1
            bar.update(count)

cur.execute("create index idx_seq on occurrence (seq)")
con.commit()
con.close()
