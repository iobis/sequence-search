import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
# import progressbar
import sqlite3


sequence_file_path = os.path.abspath("data/sequences.csv")
occurrence_file_path = os.path.abspath("data/occurrences.csv")
fasta_file_path = os.path.abspath("data/sequences.fasta")
sqlite_file_path = os.path.abspath("data/occurrence.sqlite")


# fasta file

csv.field_size_limit(100000000)

try:
    os.remove(fasta_file_path)
except FileNotFoundError:
    pass

# bar = progressbar.ProgressBar()
with open(fasta_file_path, "w") as fasta_file:
    with open(sequence_file_path) as csv_file:
        reader = csv.DictReader(csv_file)
        count = 0
        for row in reader:
            seq = str(row["sequence"]).replace("-", "")
            record = SeqRecord(
                Seq(seq),
                id=row["hash"],
                description=row["hash"]
            )
            SeqIO.write([record], fasta_file, "fasta")
            count = count + 1
            # bar.update(count)

# sqlite (occurrence IDs)

try:
    os.remove(sqlite_file_path)
except FileNotFoundError:
    pass

con = sqlite3.connect(sqlite_file_path)
cur = con.cursor()
cur.execute("create table occurrence (hash, decimallongitude real, decimallatitude real, dataset_id, phylum, class, \"order\", family, genus, scientificname, count int)")
# with progressbar.ProgressBar(max_value=progressbar.UnknownLength) as bar:
with open(occurrence_file_path) as csv_file:
    reader = csv.DictReader(csv_file)
    count = 0
    for row in reader:
        cur.execute("insert into occurrence values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", list(row.values()))
        count = count + 1
        # bar.update(count)

cur.execute("create index idx_hash on occurrence (hash)")
con.commit()
con.close()
