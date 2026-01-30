# Universal Plasmid Assembly Pipeline

A **computational pipeline for automated plasmid construction** that simulates core molecular cloning principles using genome sequences, design specifications, and biological marker libraries.

---

## Overview

This project generates a **synthetic plasmid DNA sequence** from:

- A genome or plasmid FASTA file (`.fa`)
- A design specification file (`Design.txt`)
- A library of selectable and screening marker genes

### Key Features
- **Automatic ORI detection** using GC-skew analysis
- **Modular plasmid assembly** based on user design
- **Global removal of restriction enzyme sites**
- **Tested on real plasmids** (e.g. pUC19)

---

## Biological Motivation

In real-world molecular cloning:

- Plasmids require an **origin of replication (ORI)** to propagate in host cells
- **Selection and screening** rely on:
  - Antibiotic resistance genes (AmpR, KanR, etc.)
  - Reporter genes (lacZα, GFP, etc.)
- **Restriction enzymes** define cloning boundaries

This repository provides a **computational analogue** of these processes.

---

## Repository Structure

```
├── main.py                # Pipeline entry point
├── ori_finder.py          # GC-skew based ORI detection
├── plasmid_builder.py     # Plasmid construction logic
├── restriction_sites.py   # Restriction enzyme motif database
├── markers/               # Marker gene FASTA files
│   ├── Ampicillin.fa
│   ├── Kanamycin.fa
│   └── Chloramphenicol.fa
├── markers.tab            # Supported marker registry
├── Design_pUC19.txt       # Example design file
├── pUC19.fa               # Example genome/plasmid
└── tests/
    └── test_pUC19.py      # Automated validation test
```

---

## Input Files

### Genome FASTA (`input.fa`)
Any bacterial genome or plasmid sequence.

```text
>Genome
ATGCGTAGCTAGCTAG...
```

---

### Design Specification (`Design.txt`)
Defines the components to include in the plasmid.

**Format**
```text
Label, Value
```

**Example**
```text
BamHI_site, BamHI
HindIII_site, HindIII
AmpR_gene, Ampicillin
lacZ_alpha, Blue_White_Selection
```

**Rules**
- Restriction enzymes must be defined in `restriction_sites.py`
- Marker genes must exist as FASTA files in `markers/`
- Unknown entries are skipped with warnings

---

## ORI Detection Method

ORI detection is performed using **GC-skew analysis** with a sliding window approach.

**Parameters**
- Window sizes: 150 bp, 200 bp, 250 bp, 300 bp
- Step size: 10 bp

The final ORI position is chosen as the **median start coordinate across scales**, providing:
- Robustness to noise
- No fixed ORI length assumption
- Stable consensus detection

---

## Installation & Usage

### Install dependencies
```bash
pip install biopython
```

### Run the pipeline
```bash
python main.py input.fa Design.txt
```

**Output**
- `Output.fa` — the finalized plasmid sequence

---

## Testing

Run the automated test suite:

```bash
python -m tests.test_pUC19
```

**Expected output**
```text
Test passed: EcoRI successfully removed.
```

---

## 🔬 Pipeline Workflow

1. Read genome FASTA
2. Detect ORI using GC-skew
3. Append marker genes
4. Insert restriction motifs
5. Remove all known restriction sites
6. Output final plasmid sequence

---
