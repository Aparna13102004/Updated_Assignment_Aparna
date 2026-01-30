from Bio import SeqIO
import os
from ori_finder import find_ori_multi_scale
from restriction_sites import RESTRICTION_SITES


def read_design(design_path):
    sites = list()
    antibiotic_genes = list()
    reporter_genes = list()

    with open(design_path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue

            label, value = line.split(",", 1)
            label = label.lower()
            value = value.strip()

            if "site" in label:
                sites.append(value)
            elif "gene" in label:
                antibiotic_genes.append(value)
            else:
                reporter_genes.append(value)

    return sites, antibiotic_genes, reporter_genes


def fetch_marker(marker):
    file_path = os.path.join("markers", marker + ".fa")
    if not os.path.isfile(file_path):
        print(f"[WARNING] Marker is missing, skipping: {marker}")
        return ""

    return str(SeqIO.read(file_path, "fasta").seq)


def purge_restriction_sites(sequence, enzymes):
    for enzyme in enzymes:
        motif = RESTRICTION_SITES.get(enzyme)
        if motif:
            sequence = sequence.replace(motif, "")
        else:
            print(f"[WARN] Unknown enzyme skipped: {enzyme}")
    return sequence


def build_plasmid(genome_fasta, design_file):
    # ORI detection
    ori_sequence, ori_start, ori_end = find_ori_multi_scale(genome_fasta)
    print(f"ORI region: {ori_start}-{ori_end}")


    sites, antibiotics, reporters = read_design(design_file)

    plasmid_seq = ori_sequence

    # Add antibiotic resistance genes
    for gene in antibiotics:
        plasmid_seq = plasmid_seq + fetch_marker(gene)

    # Add screening / reporter genes
    for gene in reporters:
        plasmid_seq = plasmid_seq + fetch_marker(gene)

    # Add restriction sites
    for site in sites:
        if site in RESTRICTION_SITES:
            plasmid_seq = plasmid_seq + RESTRICTION_SITES[site]
        else:
            print("[WARNING] Unknown restriction site skipping:" + site)

    # Remove all restriction motifs globally
    plasmid_seq = purge_restriction_sites(plasmid_seq, RESTRICTION_SITES.keys())

    return plasmid_seq
