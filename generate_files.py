import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import progressbar
import sqlite3


csv_file_path = os.path.abspath("data/sequences.csv")
fasta_file_path = os.path.abspath("data/sequences.fasta")
sqlite_file_path = os.path.abspath("data/occurrence.sqlite")


# fasta file

csv.field_size_limit(100000000)

try:
    os.remove(fasta_file_path)
except FileNotFoundError:
    pass

with progressbar.ProgressBar(max_value=progressbar.UnknownLength) as bar:
    with open(fasta_file_path, "w") as fasta_file:
        with open(csv_file_path) as csv_file:
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
    os.remove(sqlite_file_path)
except FileNotFoundError:
    pass

con = sqlite3.connect(sqlite_file_path)
cur = con.cursor()
cur.execute("create table occurrence (id, seq)")

with progressbar.ProgressBar(max_value=progressbar.UnknownLength) as bar:
    with open(csv_file_path) as csv_file:
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
