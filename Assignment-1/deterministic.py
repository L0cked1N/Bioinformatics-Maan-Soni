#!/usr/bin/env python3
import sys
import json

# ---------- FASTA READER ----------
def read_fasta(filepath):
    seq = []
    with open(filepath, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)

# ---------- DESIGN PARSER (UPDATED) ----------
def parse_design(design_file, restriction_sites):
    mcs_sequences = []
    marker = None

    # canonical antibiotic names we support
    known_markers = {"ampicillin", "kanamycin", "chloramphenicol"}

    with open(design_file, "r") as f:
        for line in f:
            if not line.strip():
                continue

            parts = [x.strip() for x in line.split(",")]
            if len(parts) < 2:
                raise ValueError(f"Invalid design line: {line.strip()}")

            label, value = parts[0], parts[1]
            value_clean = value.strip()
            value_lower = value_clean.lower()

            # ---- Antibiotic marker detection (ROBUST) ----
            if (
                "antibiotic" in label.lower()
                or value_lower in known_markers
            ):
                if marker is not None:
                    raise ValueError("Multiple antibiotic markers specified")
                marker = value_clean
                continue

            value_upper = value_clean.upper()

            # Literal DNA sequence
            if all(c in "ATGC" for c in value_upper):
                mcs_sequences.append(value_upper)
                continue

            # Named restriction enzyme
            if value_clean in restriction_sites:
                mcs_sequences.append(restriction_sites[value_clean])
                continue

            raise ValueError(f"Unknown restriction site or invalid DNA: {value_clean}")

    if marker is None:
        raise ValueError("No antibiotic marker specified in design.txt")

    return mcs_sequences, marker.lower()
# ---------- MAIN ----------
def main():
    if len(sys.argv) != 3:
        print("Usage: python deterministic.py input.fa design.txt")
        sys.exit(1)

    insert_fasta = sys.argv[1]
    design_file = sys.argv[2]

    # Paths based on YOUR directory structure
    backbone_path = "GenBank_Data_RSF1010/e_coli_plasmid_rsf1010.fa"

    marker_files = {
        "ampicillin": "GenBank_Data_RSF1010/Antibiotic_Resistance_Markers/Ampicillin.fa",
        "kanamycin": "GenBank_Data_RSF1010/Antibiotic_Resistance_Markers/Kanamycin.fa",
        "chloramphenicol": "GenBank_Data_RSF1010/Antibiotic_Resistance_Markers/Chloramphenicol.fa"
    }

    # Load restriction sites
    with open("restriction_sites.json", "r") as f:
        restriction_sites = json.load(f)

    # Parse design
    mcs_sequences, marker_name = parse_design(design_file, restriction_sites)

    if marker_name not in marker_files:
        raise ValueError(f"Unsupported antibiotic marker: {marker_name}")

    # Load sequences
    backbone = read_fasta(backbone_path)
    insert = read_fasta(insert_fasta)
    marker_seq = read_fasta(marker_files[marker_name])

    # ---------- ASSEMBLY (DETERMINISTIC) ----------
    plasmid = backbone

    for seq in mcs_sequences:
        plasmid += seq

    plasmid += insert
    plasmid += marker_seq

    # Write output
    with open("Output_deterministic.fa", "w") as f:
        f.write(">Constructed_Plasmid\n")
        for i in range(0, len(plasmid), 60):
            f.write(plasmid[i:i+60] + "\n")

    print("✅ Deterministic plasmid construction complete.")
    print("📄 Output written to Output_deterministic.fa")

if __name__ == "__main__":
    main()
