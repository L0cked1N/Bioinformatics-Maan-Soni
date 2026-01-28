# Universal Plasmid Construction Pipeline

## 📌 Overview

This repository implements a **universal plasmid construction pipeline** in Python that generates a complete plasmid DNA sequence from:

1. A **backbone plasmid** (RSF1010),
2. A user-provided **DNA insert** (FASTA),
3. A list of **restriction sites** and a **selected antibiotic resistance marker** defined in a design specification.

The tool comes with **two independent algorithmic solutions**:
- A **deterministic builder** (`deterministic.py`)
- An **algorithmic rule-based builder** (`algorithmic.py`)

It also includes a **comparison script** to ensure correctness and equivalence.

---

## 🧬 Biological Background

In this project:
- We use the **RSF1010** replicon as a fixed backbone.
- Antibiotic markers supported: **Ampicillin**, **Kanamycin**, **Chloramphenicol**.
- Restriction sites can be either known enzyme names (e.g., `EcoRI`) or literal DNA sequences (e.g., `GAATTC`).

---

## 📁 Repository Structure

Assignment-1/
├── deterministic.py # Deterministic plasmid builder
├── algorithmic.py # Rule-based plasmid builder
├── compare_solutions.py # Compare outputs of both builders
├── restriction_sites.json # Known restriction enzyme sites
├── input.fa # THIS WILL BE ADDED BY THE PERSON FOR TESTING AND USING THIS 
├── design.txt # THIS WILL BE ADDED BY THE PERSON FOR TESTING AND USING THIS 
├── Output.fa # Generated plasmid (after running)
├── GenBank_Data_RSF1010/
│ ├── e_coli_plasmid_rsf1010.fa # RSF1010 backbone FASTA
│ └── Antibiotic_Resistance_Markers/
│ ├── Ampicillin.fa # Ampicillin resistance marker (CDS)
│ ├── Kanamycin.fa # Kanamycin resistance marker (CDS)
│ └── Chloramphenicol.fa # Chloramphenicol resistance marker (CDS)
├── tests/
│ ├── test_plasmid_builder.py # Automatic tests
│ ├── design.txt # Test design 1
│ ├── design1.txt # Test design 2
│ ├── design2.txt # Test design 3
└── myvenv/ # Python virtual environment (optional)

---

## 🧾 Input Files
### ❗ input.fa Format
### ❗ design.txt Format

Each line is a comma-separated specification:

Label, Value

- **MCS entries**: can be named enzymes or literal DNA sequences.
  - Named enzyme → looked up in `restriction_sites.json`
  - Literal sequence → used directly (if only A/T/G/C)
- **Antibiotic marker**:
  - Must include the word “antibiotic” in the label OR a known antibiotic name as the value.
  - Only one marker may be specified.

Examples:

MCS1, EcoRI
MCS2, HindIII
Antibiotic_marker1, Ampicillin

or:

Multiple_Cloning_Site1, GAATTC
Multiple_Cloning_Site2, AAGCTT
Antibiotic_marker1, Kanamycin

( This code deals with multiple styles of input design.txt )
---

## 🧠 How It Works

### 🔹 Deterministic Builder (`deterministic.py`)

This script:

1. Reads the **backbone**, **insert**, and **antibiotic marker** sequences.
2. Parses the design to interpret restriction sites (enzyme names or literal DNA).
3. Concatenates:
backbone -> MCS sequences -> insert -> marker

4. Writes a **FASTA** file (`Output_deterministic.fa`) with the complete plasmid.

Usage:
python deterministic.py input.fa design.txt

###🔹 Algorithmic Builder (algorithmic.py)

This script takes the same inputs but:

1. Builds a component dictionary first.
2. Constructs a plan (assembly order) based on parsed design.
3. Executes the plan to generate the same plasmid sequence as the deterministic builder.

4. Writes a **FASTA** file (`Output_algo.fa`) with the complete plasmid.

Usage:
python algorithmic.py input.fa design.txt

# 📌 Assumptions & Scope
✔ Only CDS-only marker sequences (no promoters/terminators) are included.
✔ Backbone is treated as a fixed string (RSF1010 replicon).
✔ No biological simulation of restriction enzyme digestion.
✔ Outputs are simple FASTA sequences.
✔ Parser supports DNA sequences directly or enzyme names.
