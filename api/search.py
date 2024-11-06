from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pysam
import sqlite3
import tempfile
import logging


sqlite_file_path = os.path.abspath("../data/occurrence.sqlite")
logging.basicConfig(level=logging.INFO)


os.environ["PATH"] += os.pathsep + "/Users/pieter/IPOfI Dropbox/Pieter Provoost/werk/projects/OBIS pipeline/sequence-search/vsearch-2.29.1-macos-aarch64/bin"
os.environ["PATH"] += os.pathsep + "/data/vsearch-2.29.1/bin"


def search_sequences(seq, strategy="vsearch"):

    with tempfile.TemporaryDirectory() as tmpdir:

        logging.info(f"Created temp dir {tmpdir}")

        input_file_path = os.path.join(tmpdir, "input.fasta")

        # generate input file

        with open(input_file_path, "w") as fasta_file:
            record = SeqRecord(
                Seq(seq),
                id="x",
                description="input sequence"
            )
            SeqIO.write([record], fasta_file, "fasta")

        records = []

        if strategy == "bowtie2":

            sam_file_path = os.path.join(tmpdir, "output.sam")

            bowtie2_path = "/Users/pieter/IPOfI Dropbox/Pieter Provoost/werk/projects/OBIS pipeline/sequence-search/bowtie2-2.5.0-macos-arm64"
            os.environ["PATH"] += os.pathsep + bowtie2_path
            os.system(f"bowtie2 -x ../data/db/sequences -U {input_file_path} -S {sam_file_path} -f --local --very-sensitive --no-unal -k 100")

            con = sqlite3.connect(sqlite_file_path)
            con.row_factory = sqlite3.Row
            cur = con.cursor()

            samfile = pysam.AlignmentFile(sam_file_path)
            for read in samfile.fetch():

                res = cur.execute("select * from occurrence where hash = ? limit 100", (read.reference_name,))
                rows = res.fetchall()

                for row in rows:
                    record = dict(zip(row.keys(), row))
                    record["as"] = read.get_tag("AS")
                    records.append(record)

                if len(records) > 500:
                    break

            con.close()
            samfile.close()

        elif strategy == "vsearch":

            blast_file_path = os.path.join(tmpdir, "output.blast")

            os.system(f"vsearch --usearch_global {input_file_path} --db ../data/sequences.udb --id 0.9 --query_cov 0.9 --maxaccepts 100 --maxrejects 100 --maxhits 100 --blast6out {blast_file_path}")

            con = sqlite3.connect(sqlite_file_path)
            con.row_factory = sqlite3.Row
            cur = con.cursor()

            with open(blast_file_path, "r") as f:
                for line in f:
                    fields = line.strip().split("\t")
                    reference = fields[1]
                    identity = float(fields[2])

                    res = cur.execute("select * from occurrence where hash = ? limit 100", (reference,))
                    rows = res.fetchall()

                    for row in rows:
                        record = dict(zip(row.keys(), row))
                        record["as"] = identity
                        records.append(record)

                    if len(records) > 500:
                        break

        return records
