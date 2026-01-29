# plasmid_builder.py
from Bio import SeqIO
import os
from ori_finder import find_ori_multi_scale

from restriction_sites import RESTRICTION_SITES


def parse_design_file(design_file):
    mcs = []
    antibiotics = []
    screening = []

    with open(design_file) as f:
        for line in f:
            if not line.strip():
                continue

            part, name = line.strip().split(",")
            part = part.strip()
            name = name.strip()

            if "site" in part.lower():
                mcs.append(name)
            elif "gene" in part.lower():
                antibiotics.append(name)
            else:
                screening.append(name)

    return mcs, antibiotics, screening


def load_marker_sequence(marker_name):
    path = f"markers/{marker_name}.fa"
    if not os.path.exists(path):
        print(f"[WARNING] Marker not found, skipping: {marker_name}")
        return ""

    record = SeqIO.read(path, "fasta")
    return str(record.seq)


def delete_sites(seq, enzymes):
    for enz in enzymes:
        if enz not in RESTRICTION_SITES:
            print(f"[WARNING] Unknown restriction enzyme, skipping: {enz}")
            continue
        motif = RESTRICTION_SITES[enz]
        seq = seq.replace(motif, "")
    return seq


def build_plasmid(input_fasta, design_file):
    # Step 1: Find ORI in host genome
    ori_seq, start, end = find_ori_multi_scale(input_fasta)
    print(f"ORI found at: {start} - {end}")

    # Step 2: Parse design file
    mcs, antibiotics, screening = parse_design_file(design_file)

    plasmid = ori_seq

    # Step 3: Add antibiotic genes
    for gene in antibiotics:
        plasmid += load_marker_sequence(gene)

    # Step 4: Add screening genes
    for gene in screening:
        plasmid += load_marker_sequence(gene)

    # Step 5: Add MCS restriction sites
    for enzyme in mcs:
        if enzyme not in RESTRICTION_SITES:
            print(f"[WARNING] Unknown restriction enzyme, skipping: {enzyme}")
            continue
        plasmid += RESTRICTION_SITES[enzyme]

    # Step 6: Delete restriction sites globally
    plasmid = delete_sites(plasmid, RESTRICTION_SITES.keys())

    return plasmid
