from fastapi import FastAPI
from search import search

app = FastAPI()


@app.get("/")
async def root():

    seq = """
TAGTCATATGCTTGTCTCAAAGATAAGCCATGCATGTCTAAGTATAAGCGACTATACTGTGAAACTGCGA
ATGGCTCATTAAATCAGTTATGGTTTATTTGATGGTACCTTGCTACTTGGATAACCGTAGTAATTCTAGA
GCTAATACATGCAGGAGTTCCCGACTCACGGAGGGATGTATTTATTAGATAAGAAACCAAACCGGTCTCC
GGTTGCGTGCTGAGTCATAATAACTGCTCGAATCGCACGGCTCTACGCCGGCGATGGTTCATTCAAATTT
CTGCCCTATCAGCTTTCGATGGTAGGATAGAGGCCTACCATGGCGTTAACGGGTAACGGAGAATTAGGGT
TCGATTCCGGAGAGGGAGCCTGAGAAATGGCTACCACATCCAAGGAAGGCAGCAGGCGCGTAAATTGCCC
GAATCCTGACACAGGGAGGTAGTGACAAGAAATAACAATACAGGGCTATTTTAGTCTTGTAATTGGAATG
AGTACAATTTACATCTCTTCACGAGGATCAATTGGAGGGCAAGTCTGGTGCCAGCAGCCGCGGTAATTCC
AGCTCCAATAGCGTATATTAAAGTTGTTGCAGTTAAAACGCTCGTAGTCGGATTTCGGGGCGGGCCGACC
GGTCTGCCGATGGGTATGCACTGGCCGGCGCGTCCTTCCACCCGGAGACCGCGCCTACTCTTAACTGAGC
GGGCGCGGGAGACGGGTCTTTTACTTTGAAAAAATCAGAGTGTTTCAAGCAGGCAGTCGCTCTTGCATGG
ATTAGCATGGGATAATGAAATAGGACTCTGGTGCTATTTTGTTGGTTTCGAACACCGGAGTAATGATTAA
CAGGGACAGTCAGGGGCACTCGTATTCCGCCGAGAGAGGTGAAATTCTCAGACCAGCGGAAGACGAACCA
CTGCGAAAGCATTTGCCAGGGATGTTTTCACTGATCAAGAACGAAAGTTAGGGGATCGAAGACGATCAGA
TACCGTCGTAGTCTTAACCATAAACCATGCCGACTAGGGATTGGAGGATGTTCCATTTGTGACTCCTTCA
GCACCTTTCGGGAAACTAAAGTCTTTGGGTTCCGGGGGGAGTATGGTCGCAAGGCTGAAACTTAAAGGAA
TTGACGGAAGGGCACCACCAGGAGTGGAGCCTGCGGCTTAATTTGACTCAACACGGGGAAACTTACCAGG
TCCAGCACATTGTGAGGATTGACAGATTGAGAGCTCTTTCTTGATTCGATGGGTGGTGGTGCATGGCCGT
TCTTAGTTGGTGGAGTGATTTGTCTGGTTAATTCCGTTAACGAACGAGACCGCAGCCTGCTAAATAGCGA
CGCGAACCCTCCGTTCGCTGGAGCTTCTTAGAGGGACAACTTGTCTTCAACAAGTGGAAGTTCGCGGCAA
TAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGATGCACTCAACGAGTCTA
TCACCTTGACCGAGAGGTCCGGGTAATCTTTTGAAATTGCATCGTGATGGGGATAGATTATTGCAACTAT
TAATCTTCAACGAGGAATTCCTAGTAAGCGTGTGTCATCAGCGCACGTTGATTACGTCCCTGCCCTTTGT
ACACACCGCCCGTCGCTCCTACCGATTGAATGATCCGGTGAGGCCCCCGGACTGCGGCGCCGCAGCTGGT
TCTCCAGCCGCGACGCCGCGGGAAGCTGTCCGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAAGG
TTTCC
    """.replace("\n", "")
    results = search(seq)

    return {"results": results}