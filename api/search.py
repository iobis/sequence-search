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


def search_sequences(seq):

    with tempfile.TemporaryDirectory() as tmpdir:

        logging.info(f"Created temp dir {tmpdir}")

        input_file_path = os.path.join(tmpdir, "input.fasta")
        sam_file_path = os.path.join(tmpdir, "output.sam")

        # generate input file

        with open(input_file_path, "w") as fasta_file:
            record = SeqRecord(
                Seq(seq),
                id="x",
                description="input sequence"
            )
            SeqIO.write([record], fasta_file, "fasta")

        # run bowtie2

        os.system(f"bowtie2 -x ../data/db/sequences -U {input_file_path} -S {sam_file_path} -f --local --very-sensitive --no-unal -k 100")

        # sam results

        records = []

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

        return records
