from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os


seq = """
TTAAAGTTGTTGTGGTTAAAAAGCTCGTAGTTGGATCTCAACAGACCTGGATCGGTTTACTTTGTGTAATACTGATCCAGAGCTGTGCTTTTGCCGGAGGTTTAGGGGTGCTCTTAACCGAGTGTCTCTGAGTGCCGGCAGGTTTACTTTGAAAAAATTAGAGTGTTTAAAGCTGTTATGCCTGAATATTTGTGCATGGAATAATAGAATAGGATGTTGATCTTATTTTGTTGGTTTTCGAGACTTGACATAATGATTAATAGGGACAGTCGGCGGCATTTGTATTCAAACGACAGAGGTGAAATTCTTGGACCGTTTGAAGACAAACTACTGCGA
""".replace("\n", "")

try:
    os.remove("input.fastq")
    os.remove("input.fasta")
except FileNotFoundError:
    pass

with open("input.fastq", "w") as fastq_file:
    record = SeqRecord(
        Seq(seq),
        id="x",
        description="input sequence"
    )
    record.letter_annotations["phred_quality"] = [40] * len(record.seq)
    SeqIO.write([record], fastq_file, "fastq")

with open("input.fasta", "w") as fasta_file:
    record = SeqRecord(
        Seq(seq),
        id="x",
        description="input sequence"
    )
    SeqIO.write([record], fasta_file, "fasta")
