#!/usr/bin/env python3
import subprocess
import sys

DET_CMD = ["python", "deterministic.py", "input.fa", "design.txt"]
ALG_CMD = ["python", "algorithmic.py", "input.fa", "design.txt"]

def read_output_deterministic():
    with open("Output_deterministic.fa") as f:
        return "".join(
            line.strip() for line in f if not line.startswith(">")
        )
def read_output_algo():
    with open("Output_algo.fa") as f:
        return "".join(
            line.strip() for line in f if not line.startswith(">")
        )
def main():
    print("▶ Running deterministic solution...")
    subprocess.run(DET_CMD, check=True)
    det_seq = read_output_deterministic()

    print("▶ Running algorithmic solution...")
    subprocess.run(ALG_CMD, check=True)
    alg_seq = read_output_algo()

    print("\n--- Comparison Result ---")

    if det_seq == alg_seq:
        print("✅ PASS: Deterministic and Algorithmic outputs are IDENTICAL")
        print(f"Plasmid length: {len(det_seq)} bp")
    else:
        print("❌ FAIL: Outputs are DIFFERENT")
        print(f"Deterministic length: {len(det_seq)} bp")
        print(f"Algorithmic length:   {len(alg_seq)} bp")

        # Optional: show first difference
        for i, (a, b) in enumerate(zip(det_seq, alg_seq)):
            if a != b:
                print(f"First difference at position {i}")
                print(f"Deterministic: {a}")
                print(f"Algorithmic:   {b}")
                break

        sys.exit(1)

if __name__ == "__main__":
    main()
