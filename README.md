# Universal Plasmid Assembly Pipeline

A **computational pipeline for automated plasmid construction** that simulates core molecular cloning principles using **motif-based ORI discovery**, user-defined plasmid designs, and biological marker libraries.

---

## Overview

This project generates a **synthetic plasmid DNA sequence** from:

- A genome or plasmid FASTA file (`.fa`)
- A design specification file (`Design.txt`)
- A library of selectable and screening marker genes

Unlike naive implementations, this pipeline **computationally discovers the origin of replication (ORI)** from the input sequence before constructing a plasmid backbone.

---

## Key Features

- **MEME-like ORI detection** using:
  - Position Weight Matrices (PWM)
  - Log-likelihood scoring
  - Background nucleotide modeling
- **Variable motif width discovery** (k = 6–15)
- **Design-driven plasmid assembly**
- **Biologically correct restriction site handling**
  - Restriction sites removed from functional regions
  - Multiple Cloning Site (MCS) appended last
- **Graceful handling of missing markers**
- **Validated on real plasmids (pUC19)**

---

## Biological Motivation

In real molecular cloning:

- Plasmids must contain an **origin of replication (ORI)** compatible with the host
- ORIs are defined by **sequence motifs**, not GC skew alone
- **Selection and screening** are achieved using:
  - Antibiotic resistance genes (AmpR, KanR, etc.)
  - Reporter genes (lacZα, GFP, etc.)
- **Restriction enzymes** define cloning boundaries and insertion sites

This repository provides a **computational analogue** of these biological principles.

---

## Repository Structure

```
Updated_Assignment_Aparna/
├── plasmid_builder.py      # CLI entry point + plasmid construction
├── ori_finder.py           # MEME-like ORI discovery (PWM + likelihood)
├── restriction_sites.py    # Restriction enzyme motif database
├── markers/                # Marker gene FASTA files
│   ├── Ampicillin.fa
│   ├── Kanamycin.fa
│   └── ...
├── Design_pUC19.txt        # Example design file
├── pUC19.fa                # Test plasmid sequence
├── Output.fa               # Generated plasmid (after execution)
└── tests/
    └── test_pUC19.py       # Automated validation test
```

---

## Input Files

### 1. Genome / Plasmid FASTA (`input.fa`)

Any bacterial genome or plasmid sequence.

```text
>Genome
ATGCGTAGCTAGCTAG...
```

---

### 2. Design Specification (`Design.txt`)

Defines plasmid components to be incorporated.

#### Format
```text
Label, Value
```

#### Example
```text
BamHI_site, BamHI
HindIII_site, HindIII
AmpR_gene, Ampicillin
lacZ_alpha, Blue_White_Selection
```

#### Rules
- Restriction enzymes must exist in `restriction_sites.py`
- Marker genes must exist as FASTA files in `markers/`
- Missing or unknown entries are skipped with warnings

---

## ORI Detection Method (MEME-like)

ORI discovery is performed using a **motif enrichment and probabilistic scoring strategy**, inspired by the MEME motif discovery algorithm.

### Key Characteristics

- Uses **Position Weight Matrices (PWM)** instead of raw k-mer counting
- Scores motifs using **log-likelihood ratios**
- Models **background nucleotide distribution**
- Allows **variable motif lengths (k = 6–15)**
- Prefers **statistically stable longer motifs** when supported by data
- Uses a **fixed sliding window (500 bp)** for ORI localization

This approach is **conceptually aligned with MEME**, while remaining computationally tractable for coursework.

---

## Installation & Usage

### Install dependencies
```bash
pip install biopython
```

### Run the pipeline
```bash
python plasmid_builder.py input.fa Design.txt
```

### Output
- `Output.fa` — the finalized plasmid DNA sequence
- Console output — MEME-style ORI discovery summary, including:
  - ORI coordinates
  - ORI sequence
  - Optimal motif width (k)
  - Mismatch policy

---

## Example Test Case (pUC19)

The pipeline has been validated using:

- `pUC19.fa` as input
- `Design_pUC19.txt` as design specification

### Expected behavior
- ORI detected automatically
- EcoRI restriction site removed from functional regions
- Final plasmid written to `Output.fa`

Verification:
```bash
grep GAATTC Output.fa
```
(No output confirms successful EcoRI removal.)

---

## Pipeline Workflow

1. Read genome / plasmid FASTA
2. Estimate background nucleotide frequencies
3. Scan sequence using MEME-like PWM scoring
4. Identify optimal ORI region and motif width
5. Use ORI as replication backbone
6. Append marker genes
7. Remove restriction sites from functional DNA
8. Append Multiple Cloning Site (MCS)
9. Output final plasmid sequence

---

## Notes on Design Choices

- This implementation is **MEME-inspired**, not a full EM-based MEME reimplementation
- Sliding windows are retained for ORI localization
- PWM scoring reduces bias toward short exact motifs
- Parameter choices prioritize **biological plausibility over overfitting**
