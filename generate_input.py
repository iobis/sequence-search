from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os


seq = """
TTAAAGTTGTTGTGGTTAAAAAGCTCGTAGTTGGATCTCAACAGACCTGGATCGGTTTACTTTGTGTAATACTGATCCAGAGCTGTGCTTTTGCCGGAGGTTTAGGGGTGCTCTTAACCGAGTGTCTCTGAGTGCCGGCAGGTTTACTTTGAAAAAATTAGAGTGTTTAAAGCTGTTATGCCTGAATATTTGTGCATGGAATAATAGAATAGGATGTTGATCTTATTTTGTTGGTTTTCGAGACTTGACATAATGATTAATAGGGACAGTCGGCGGCATTTGTATTCAAACGACAGAGGTGAAATTCTTGGACCGTTTGAAGACAAACTACTGCGA
""".replace("\n", "")

try:
    os.remove("input.fasta")
except FileNotFoundError:
    pass

with open("input.fasta", "w") as fasta_file:
    record = SeqRecord(
        Seq(seq),
        id="x",
        description="input sequence"
    )
    SeqIO.write([record], fasta_file, "fasta")
