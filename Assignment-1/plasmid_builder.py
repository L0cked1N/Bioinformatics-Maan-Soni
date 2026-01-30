# plasmid_builder.py
from Bio import SeqIO
import os
from ori_finder import find_ori_meme_style
from restriction_sites import RESTRICTION_SITES

def parse_design_file(design_file):
    mcs, antibiotics, screening = [], [], []
    with open(design_file) as f:
        for line in f:
            if not line.strip():
                continue
            part, name = line.strip().split(",")
            part = part.strip().lower()
            name = name.strip()      

            if "site" in part:
                mcs.append(name)
            elif "gene" in part:
                antibiotics.append(name)
            else:
                screening.append(name)
    return mcs, antibiotics, screening

def load_marker_sequence(marker_name):
    path = f"markers/{marker_name}.fa"
    if not os.path.exists(path):
        print(f"[WARNING] Marker not found: {marker_name}")
        return ""
    record = SeqIO.read(path, "fasta")
    return str(record.seq).upper()


def delete_sites(seq, enzymes):
    for enz in enzymes:
        if enz in RESTRICTION_SITES:
            seq = seq.replace(RESTRICTION_SITES[enz], "")
    return seq


def build_plasmid(input_fasta, design_file):
    # ---------- ORI DISCOVERY ----------
    ori = find_ori_meme_style(input_fasta)

    ori_seq = ori["ori_sequence"]
    start, end = ori["ori_region"]

    print("\n[ORI DISCOVERY REPORT]")
    print(f"Position        : {ori['ori_position']}")
    print(f"Region          : {start}-{end}")
    print(f"Motif length    : {ori['best_k']}")
    print(f"Log-likelihood  : {round(ori['log_likelihood'],2)}")
    print(f"Motifs          : {ori['motifs']}\n")

    # ---------- DESIGN FILE ----------
    mcs, antibiotics, screening = parse_design_file(design_file)

    # ---------- BUILD BACKBONE ----------
    plasmid = ori_seq

    for gene in antibiotics:
        plasmid += load_marker_sequence(gene)

    for gene in screening:
        plasmid += load_marker_sequence(gene)

    # ---------- SANITIZE BACKBONE ----------
    # remove restriction sites ONLY from functional DNA
    plasmid = delete_sites(plasmid, RESTRICTION_SITES.keys())

    # ---------- ADD MCS LAST ----------
    # these must NEVER be deleted
    for enzyme in mcs:
        if enzyme not in RESTRICTION_SITES:
            print(f"[WARNING] Unknown restriction enzyme: {enzyme}")
            continue
        plasmid += RESTRICTION_SITES[enzyme]

    return plasmid
