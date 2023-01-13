from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pysam
import sqlite3
import requests


input_file_path = os.path.abspath("../data/input.fasta")
sqlite_file_path = os.path.abspath("../data/occurrence.sqlite")
sam_file_path = os.path.abspath("../data/output.sam")


def search(seq):

    # cleanup

    cleanup()

    # generate input file

    with open(input_file_path, "w") as fasta_file:
        record = SeqRecord(
            Seq(seq),
            id="x",
            description="input sequence"
        )
        SeqIO.write([record], fasta_file, "fasta")

    # run bowtie2

    os.system("../bowtie2-2.5.0-macos-arm64/bowtie2 -x ../data/db/sequences -U ../data/input.fasta -S ../data/output.sam -f --local --very-sensitive --no-unal -k 100")

    # sam results

    records = []

    con = sqlite3.connect(sqlite_file_path)
    cur = con.cursor()

    samfile = pysam.AlignmentFile(sam_file_path)
    for read in samfile.fetch():

        res = cur.execute("select id from occurrence where seq = ?", (read.reference_name,))
        rows = res.fetchall()
        for row in rows:
            record = {
                "as": read.get_tag("AS"),
                "reference_name": read.reference_name,
                "query_alignment_start": read.query_alignment_start,
                "query_alignment_end": read.query_alignment_end,
                "query_alignment_length": read.query_alignment_length,
                "reference_start": read.reference_start,
                "reference_end": read.reference_end,
                "reference_length": read.reference_length,
                "occurrence_id": row[0]
            }
            records.append(record)

    con.close()
    samfile.close()

    # add OBIS occurrences

    occurrence_ids = ",".join([record["occurrence_id"] for record in records])
    url = f"https://api.obis.org/occurrence?size=1000&id={occurrence_ids}"
    res = requests.get(url)
    occurrences = res.json()["results"]

    occurrence_map = {occurrence["id"]: {
        "scientificName": occurrence["scientificName"] if "scientificName" in occurrence else None,
        "decimalLongitude": occurrence["decimalLongitude"] if "decimalLongitude" in occurrence else None,
        "decimalLatitude": occurrence["decimalLatitude"] if "decimalLatitude" in occurrence else None,
        "eventDate": occurrence["eventDate"] if "eventDate" in occurrence else None,
        "date_year": occurrence["date_year"] if "date_year" in occurrence else None,
        "dataset_id": occurrence["dataset_id"] if "dataset_id" in occurrence else None,
        "phylum": occurrence["phylum"] if "phylum" in occurrence else None,
        "class": occurrence["class"] if "class" in occurrence else None,
        "order": occurrence["order"] if "order" in occurrence else None,
        "family": occurrence["family"] if "family" in occurrence else None,
        "genus": occurrence["genus"] if "genus" in occurrence else None
    } for occurrence in occurrences}

    for record in records:
        if record["occurrence_id"] in occurrence_map:
            record.update(occurrence_map[record["occurrence_id"]])

    # cleanup

    cleanup()

    return records


def cleanup():
    try:
        os.remove(input_file_path)
    except FileNotFoundError:
        pass
    try:
        os.remove(sam_file_path)
    except FileNotFoundError:
        pass
