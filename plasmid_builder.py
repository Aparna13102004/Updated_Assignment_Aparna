"""
PLASMID BUILDER MODULE (Final Assignment Version)

Usage:
------
python plasmid_builder.py <input_fasta.fa> <design.txt>

Workflow:
---------
1. Detect ORI from unknown organism using MEME-inspired motif enrichment
2. Use detected ORI as replication backbone
3. Parse design file to determine:
   - Multiple cloning sites (restriction enzymes)
   - Antibiotic resistance markers
   - Screening / reporter genes
4. Remove restriction sites from functional DNA
5. Append MCS sites at the end (never deleted)
6. Output plasmid sequence as Output.fa
"""

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ori_finder import find_ori_multi_scale
from restriction_sites import RESTRICTION_SITES


# ---------------------------------------------------------
# Parse Design.txt (assignment format)
# ---------------------------------------------------------
def parse_design_file(design_file):
    """
    Expected Design.txt format:

    Multiple_Cloning_Site1, RestrictionEnzyme1
    Multiple_Cloning_Site2, RestrictionEnzyme2
    ...
    Antibiotic_marker1, MarkerName1
    Antibiotic_marker2, MarkerName2
    """

    mcs = []
    antibiotics = []
    screening = []

    with open(design_file) as f:
        for line in f:
            if not line.strip():
                continue

            label, value = line.strip().split(",", 1)
            label = label.lower()
            value = value.strip()

            if "site" in label:
                mcs.append(value)
            elif "antibiotic" in label or "marker" in label or "gene" in label:
                antibiotics.append(value)
            else:
                screening.append(value)

    return mcs, antibiotics, screening


# ---------------------------------------------------------
# Load marker gene sequence
# ---------------------------------------------------------
def load_marker_sequence(marker_name):
    """
    Loads marker sequence from markers/<marker_name>.fa
    Gracefully skips missing markers as required by assignment.
    """
    path = os.path.join("markers", marker_name + ".fa")

    if not os.path.exists(path):
        print(f"[WARNING] Marker not found, skipping: {marker_name}")
        return ""

    record = SeqIO.read(path, "fasta")
    return str(record.seq).upper()


# ---------------------------------------------------------
# Remove restriction sites from functional DNA
# ---------------------------------------------------------
def delete_restriction_sites(sequence, enzymes):
    """
    Removes restriction enzyme recognition sites
    from ORI and gene regions.
    """
    for enzyme in enzymes:
        motif = RESTRICTION_SITES.get(enzyme)
        if motif:
            sequence = sequence.replace(motif, "")
    return sequence


# ---------------------------------------------------------
# MEME-style ORI summary printer
# ---------------------------------------------------------
def print_meme_style_summary(ori_seq, start, end, best_k):
    """
    Prints MEME-style ORI discovery summary using
    true detected parameters.
    """

    print("\n========== ORI DISCOVERY SUMMARY (MEME-style) ==========")
    print(f"ORI length             : {end - start} bp")
    print(f"ORI coordinates        : {start} - {end}")
    print(f"Window size            : 500 bp")
    print(f"Optimal k-mer length   : {best_k}")

    print("\nMismatch policy:")
    if best_k <= 5:
        print("  d = 0-1  (short, strict motifs)")
    elif best_k <= 9:
        print("  d = 0-2  (moderate degeneracy)")
    else:
        print("  d = 0-3  (long, degenerate motifs)")

    print("\nPredicted ORI sequence:")
    print(ori_seq)

    print("========================================================\n")


# ---------------------------------------------------------
# Core plasmid construction logic
# ---------------------------------------------------------
def build_plasmid(input_fasta, design_file):
    """
    Builds plasmid using:
    ORI -> antibiotic markers -> screening genes
    -> restriction site sanitization -> MCS addition
    """

    # ---------- ORI DISCOVERY ----------
    ori_seq, start, end, best_k = find_ori_multi_scale(input_fasta)

    print_meme_style_summary(
        ori_seq=ori_seq,
        start=start,
        end=end,
        best_k=best_k
    )

    # ---------- DESIGN FILE ----------
    mcs, antibiotics, screening = parse_design_file(design_file)

    # ---------- BUILD BACKBONE ----------
    plasmid = ori_seq

    for gene in antibiotics:
        plasmid += load_marker_sequence(gene)

    for gene in screening:
        plasmid += load_marker_sequence(gene)

    # ---------- SANITIZE FUNCTIONAL DNA ----------
    plasmid = delete_restriction_sites(plasmid, RESTRICTION_SITES.keys())

    # ---------- ADD MCS LAST ----------
    for enzyme in mcs:
        if enzyme not in RESTRICTION_SITES:
            print(f"[WARNING] Unknown restriction enzyme skipped: {enzyme}")
            continue
        plasmid += RESTRICTION_SITES[enzyme]

    return plasmid


# ---------------------------------------------------------
# CLI ENTRY POINT
# ---------------------------------------------------------
if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python plasmid_builder.py <input.fa> <design.txt>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    design_file = sys.argv[2]

    plasmid_sequence = build_plasmid(input_fasta, design_file)

    record = SeqRecord(
        Seq(plasmid_sequence),
        id="Output_Plasmid",
        description="Plasmid constructed using detected ORI and design file"
    )

    SeqIO.write(record, "Output.fa", "fasta")
    print("[SUCCESS] Plasmid sequence written to Output.fa")
