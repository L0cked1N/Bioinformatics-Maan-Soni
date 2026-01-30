# tests/test_pUC19.py

from plasmid_builder import build_plasmid

def test_pUC19_restriction_site_removal():
    plasmid = build_plasmid("pUC19.fa", "Design_pUC19.txt")

    # EcoRI site must be removed
    assert "GAATTC" not in plasmid, "EcoRI site still present!"

    # Basic sanity checks
    assert plasmid.count("G") + plasmid.count("C") > 0, "Invalid DNA sequence"

    print("Test passed: EcoRI successfully removed.")


if __name__ == "__main__":
    test_pUC19_restriction_site_removal()
