#!/usr/bin/env python3
"""
Component 2: Real-Genome Validation
====================================
Reproduces Table 3 of:
  "Minimum Sequencing Depth Guarantees for Reference-Based Detection
   of Divergent Pathogens in Clinical Metagenomic Sequencing"

Each scenario applies a biologically motivated divergence mu to a real
reference genome and measures the empirical detection rate, which is then
compared against D_theory = 1 - (1 - p*q)^T from Theorem 3.1.

Note on the two cross-family scenarios: the divergence is applied to the
reference listed in the `reference` field below, not to a second downloaded
genome. The SARS-CoV-1 row models a SARS-CoV-1-like pathogen (~20% divergent
from SARS-CoV-2) by mutating the Wuhan reference at mu = 0.20, and the
cross-subtype row models an H1N1-like HA (~30% divergent from H3N2 HA) by
mutating the H3N2 HA reference at mu = 0.30. data/h1n1_ha.fasta is included
as the source of that 30% estimate; it is not itself a simulation reference,
which is why the cross-subtype row reports N = 1,737 (H3N2 HA) and not 1,752.

Usage (from the repository root):
    python code/real_genomes.py
    python code/real_genomes.py --trials 400 --seed 7
"""

import argparse
import csv
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from simulator import (  # noqa: E402  (path set above)
    load_fasta,
    q_recurrence,
    detection_prob_theory,
    simulate_detection,
    mu_c_numerical,
)

# label, reference FASTA, divergence
SCENARIOS = [
    ("SARS-CoV-2 Wuhan→Wuhan",   "sars_cov2_wuhan.fasta", 0.0000),
    ("SARS-CoV-2 Wuhan→Omicron", "sars_cov2_wuhan.fasta", 0.0017),
    ("SARS-CoV-2 Wuhan→SARS-1",  "sars_cov2_wuhan.fasta", 0.2000),
    ("H3N2 HA: same season",          "h3n2_ha.fasta",         0.0000),
    ("H3N2 HA: 1-year drift",         "h3n2_ha.fasta",         0.0100),
    ("H3N2 HA: 5-year drift",         "h3n2_ha.fasta",         0.0500),
    ("H3N2 vs H1N1 HA",               "h3n2_ha.fasta",         0.3000),
]


def main():
    parser = argparse.ArgumentParser(
        description="Real-genome validation — reproduces paper Table 3",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--datadir', default='data')
    parser.add_argument('--p',      type=float, default=1e-4)
    parser.add_argument('--w',      type=int,   default=24)
    parser.add_argument('--L',      type=int,   default=150)
    parser.add_argument('--delta',  type=float, default=0.01)
    parser.add_argument('--depth',  type=int,   default=100_000)
    parser.add_argument('--trials', type=int,   default=2000)
    parser.add_argument('--seed',   type=int,   default=42)
    parser.add_argument('--out',    default='data/results_real_genomes.csv')
    args = parser.parse_args()

    random.seed(args.seed)

    mu_c = mu_c_numerical(args.w, args.L, args.p, args.delta, args.depth)
    HR = '─' * 92
    print(f"\n{HR}")
    print("  REAL-GENOME VALIDATION (paper Table 3)")
    print(f"  p={args.p:.1e}  w={args.w}  L={args.L}  T={args.depth:,}  "
          f"δ={args.delta}  trials={args.trials}  seed={args.seed}")
    print(f"  μ_c = {mu_c:.4f} (at T_max = T)")
    print(HR)
    print(f"  {'Scenario':<30} {'N':>7}  {'μ':>6}  {'q_exact':>9}  "
          f"{'D_theory':>9}  {'D_emp':>7}  {'SE':>7}  {'|Δ|':>7}  {'Tier':>4}")
    print(HR)

    # Cache genomes so each FASTA is parsed once.
    genomes = {}
    rows = []
    for label, fasta, mu in SCENARIOS:
        if fasta not in genomes:
            genomes[fasta] = load_fasta(os.path.join(args.datadir, fasta))
        ref = genomes[fasta]

        qe = q_recurrence(args.w, args.L, mu)
        d_th = detection_prob_theory(args.p, mu, args.w, args.L, args.depth)
        d_em = simulate_detection(ref, mu, args.p, args.w, args.L,
                                  args.depth, args.trials)
        gap = abs(d_em - d_th)
        # Binomial standard error of the Monte Carlo estimate itself, so the
        # reader can judge which gaps are sampling noise and which are the
        # shared-genome (Jensen) effect of Section 5.4.
        se = math.sqrt(d_em * (1.0 - d_em) / args.trials)
        tier = '1' if mu == 0.0 else ('2' if mu <= mu_c else '3')

        print(f"  {label:<30} {len(ref):>7,}  {mu:>6.4f}  {qe:>9.5f}  "
              f"{d_th:>9.4f}  {d_em:>7.4f}  {se:>7.4f}  {gap:>7.4f}  {tier:>4}")
        rows.append(dict(label=label, N=len(ref), mu=mu, q_exact=qe,
                         det_theory=d_th, det_empirical=d_em, se=se, gap=gap,
                         tier=tier, reference=fasta, trials=args.trials,
                         seed=args.seed))
    print(HR)
    print(f"  Max |Δ| across real genomes: {max(r['gap'] for r in rows):.4f}")

    # encoding is explicit: the scenario labels contain "→", which the default
    # Windows cp1252 codec cannot encode.
    with open(args.out, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            'label', 'N', 'mu', 'q_exact', 'det_theory', 'det_empirical',
            'se', 'gap', 'tier', 'reference', 'trials', 'seed'
        ])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nResults → {args.out}")
    print("Done.")


if __name__ == '__main__':
    main()
