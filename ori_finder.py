"""
ORI FINDER MODULE (MEME-like PWM + Log-Likelihood)

This implementation approximates MEME by:
- Building position weight matrices (PWM)
- Scoring motifs using log-likelihood ratios
- Allowing variable motif widths (k)
- Using background nucleotide frequencies
"""

from Bio import SeqIO
import math
from collections import defaultdict

NUCLEOTIDES = ["A", "C", "G", "T"]


# ---------------------------------------------------------
# Compute background nucleotide distribution
# ---------------------------------------------------------
def background_freq(sequence):
    counts = defaultdict(int)
    for nt in sequence:
        counts[nt] += 1
    total = sum(counts.values())
    return {nt: counts[nt] / total for nt in NUCLEOTIDES}


# ---------------------------------------------------------
# Build PWM for a set of aligned motifs
# ---------------------------------------------------------
def build_pwm(motifs, pseudocount=1):
    k = len(motifs[0])
    pwm = [{nt: pseudocount for nt in NUCLEOTIDES} for _ in range(k)]

    for motif in motifs:
        for i, nt in enumerate(motif):
            pwm[i][nt] += 1

    # Normalize to probabilities
    for i in range(k):
        total = sum(pwm[i].values())
        for nt in NUCLEOTIDES:
            pwm[i][nt] /= total

    return pwm


# ---------------------------------------------------------
# Compute log-likelihood of motifs under PWM
# ---------------------------------------------------------
def pwm_log_likelihood(motifs, pwm, bg):
    score = 0.0
    for motif in motifs:
        for i, nt in enumerate(motif):
            score += math.log(pwm[i][nt] / bg[nt])
    return score


# ---------------------------------------------------------
# MEME-like motif scoring for a window
# ---------------------------------------------------------
def meme_style_score(window, k, bg):
    motifs = [window[i:i+k] for i in range(len(window) - k + 1)]

    if len(motifs) < 5:
        return -float("inf")

    pwm = build_pwm(motifs)
    return pwm_log_likelihood(motifs, pwm, bg)


# ---------------------------------------------------------
# ORI scanning using MEME-like logic
# ---------------------------------------------------------
def scan_ori_region(
    fasta_path,
    window_size=500,
    step_size=25,
    k_range=range(6, 16)
):
    record = SeqIO.read(fasta_path, "fasta")
    sequence = str(record.seq).upper()

    bg = background_freq(sequence)

    best_score = -float("inf")
    best_start = 0
    best_k = None

    for start in range(0, len(sequence) - window_size, step_size):
        window = sequence[start:start + window_size]

        for k in k_range:
            score = meme_style_score(window, k, bg)

            if score > best_score:
                best_score = score
                best_start = start
                best_k = k

    return best_start, best_start + window_size, best_k


# ---------------------------------------------------------
# Public API
# ---------------------------------------------------------
def find_ori_multi_scale(fasta_path):
    """
    Returns:
    - ORI sequence
    - start
    - end
    - best motif width (k)
    """
    start, end, best_k = scan_ori_region(fasta_path)

    record = SeqIO.read(fasta_path, "fasta")
    sequence = str(record.seq).upper()

    return sequence[start:end], start, end, best_k
