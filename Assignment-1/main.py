# main.py
import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from plasmid_builder import build_plasmid

if len(sys.argv) != 3:
    print("Usage: python main.py Input.fa Design.txt")
    sys.exit(1)

input_fasta = sys.argv[1]
design_file = sys.argv[2]

plasmid_seq = build_plasmid(input_fasta, design_file)

record = SeqRecord(
    Seq(plasmid_seq),
    id="Output_Plasmid",
    description="Plasmid constructed using host ORI and design file"
)

SeqIO.write(record, "Output.fa", "fasta")
print("Plasmid written to Output.fa")
