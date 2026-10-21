#!/usr/bin/env python3
"""
Component 1: Metagenomic Read Simulator
========================================
Validates Theorems 3.1, 3.3 and 3.4 from:
  "Minimum Sequencing Depth Guarantees for Reference-Based Detection
   of Divergent Pathogens in Clinical Metagenomic Sequencing"

Each trial generates a fresh mutated genome so the empirical estimate
reflects the expectation over all possible diverged sequences, not one
specific realisation. Within a trial, however, every read comes from the
same mutated genome; see Section 5.4 and code/shared_genome_expectation.py
for the resulting bias, which is of known sign and leaves T* conservative.

Usage (from the repository root):
    python code/simulator.py                       # reproduces paper Tables 1-2
    python code/simulator.py --genome my.fasta     # real FASTA
    python code/simulator.py --p 1e-4 --depth 200000 --trials 200 --seed 42
"""

import argparse
import csv
import math
import random
import sys

# ─── Sequence utilities ───────────────────────────────────────────────────────

BASES = ['A', 'T', 'G', 'C']


def random_base(exclude=None):
    choices = [b for b in BASES if b != exclude]
    return random.choice(choices)


def random_sequence(length):
    return ''.join(random.choices(BASES, k=length))


def load_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith('>'):
                seq.append(line.upper().replace(' ', ''))
    return ''.join(seq)


# ─── Mutation model ───────────────────────────────────────────────────────────

def mutate_sequence(seq, mu):
    """Independent substitution mutations at per-base rate mu."""
    if mu == 0.0:
        return seq
    result = []
    for base in seq:
        if random.random() < mu:
            result.append(random_base(exclude=base))
        else:
            result.append(base)
    return ''.join(result)


# ─── Read simulation ──────────────────────────────────────────────────────────

def sample_read(genome, L):
    pos = random.randint(0, len(genome) - L)
    return genome[pos:pos + L], pos


def is_informative(read, reference_region, w):
    """True if read has at least one w-mer exactly matching the reference."""
    for i in range(len(read) - w + 1):
        if read[i:i + w] == reference_region[i:i + w]:
            return True
    return False


# ─── Theoretical bounds ───────────────────────────────────────────────────────

def q_recurrence(w, L, mu):
    """
    Exact informative-read probability (Definition 2.3, Eqs. 2-4).
    P[n] = P(no run of w consecutive matches in n Bernoulli(1-mu) trials).
    Uses linear recurrence of order w; O(wL) time.

    The boundary condition P[n] = 1 for n < w is essential. It is applied here
    via the explicit all-successes term below, which makes the n < w rows sum
    to 1; dropping it gives P[1] = mu and corrupts every q (see Remark 2.5).
    """
    p = 1.0 - mu
    if L < w:
        return 0.0
    P = [0.0] * (L + 1)
    P[0] = 1.0
    for n in range(1, L + 1):
        total = 0.0
        for j in range(min(n, w)):
            idx = n - j - 1
            total += (p ** j) * mu * P[idx]
        if n < w:
            total += p ** n   # all-successes case, safe since n < w
        P[n] = total
    return 1.0 - P[L]


def q_poisson(w, L, mu):
    """
    Poisson approximation (Eq. 5). For comparison only -- it overestimates q
    because overlapping matching windows clump: a run of m >= w matches yields
    m-w+1 matching windows, so lam overcounts independent chances by ~1/mu.
    A small lam does not rescue it (at mu=0.30, lam ~ 0.024 yet the ratio to
    the exact value is still 3.2). Use q_recurrence for any real planning.
    """
    lam = (L - w + 1) * ((1.0 - mu) ** w)
    return 1.0 - math.exp(-lam)


def T_star(p, mu, w, L, delta, use_exact=True):
    """Minimum depth for delta-detection (Theorem 3.1, Eq. 6)."""
    q = q_recurrence(w, L, mu) if use_exact else q_poisson(w, L, mu)
    if q <= 0.0:
        return float('inf')
    pq = p * q
    if pq <= 0.0:
        return float('inf')
    log_term = math.log(1.0 - pq)
    if log_term >= 0:
        return 1
    return math.ceil(math.log(delta) / log_term)


def detection_prob_theory(p, mu, w, L, T):
    """P(detection) = 1 - (1-p*q)^T from Theorem 3.1, Eq. 6."""
    q = q_recurrence(w, L, mu)
    pq = p * q
    if pq <= 0.0:
        return 0.0
    return 1.0 - (1.0 - pq) ** T


def mu_c_numerical(w, L, p, delta, T_max, tol=1e-7):
    """Bisection for critical divergence mu_c (Theorem 3.3, Eq. 7)."""
    if T_max <= 0 or T_max == float('inf'):
        return 0.0
    target = math.log(1.0 / delta) / (p * T_max)
    if target >= 1.0:
        return 0.0
    if q_recurrence(w, L, 0.0) < target:
        return 0.0
    lo, hi = 0.0, 1.0 - 1e-9
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if q_recurrence(w, L, mid) > target:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return (lo + hi) / 2.0


def w_star(Q, n_avg, alpha):
    """Optimal window length (Theorem 3.4, Eq. 9)."""
    return math.ceil(math.log(Q * n_avg / alpha) / math.log(4))


# ─── Core simulation engine ───────────────────────────────────────────────────

def simulate_detection(reference, mu, p, w, L, total_reads, n_trials):
    """
    Run n_trials independent experiments.
    KEY FIX: each trial generates a fresh mutated pathogen genome,
    so the empirical rate averages over the distribution of all
    possible mutated sequences — not one specific realisation.
    """
    N = len(reference)
    detected = 0

    for _ in range(n_trials):
        # Fresh mutation draw per trial
        pathogen = mutate_sequence(reference, mu)
        found = False
        for _ in range(total_reads):
            if random.random() < p:
                read, pos = sample_read(pathogen, L)
                ref_region = reference[pos:pos + L]
                if is_informative(read, ref_region, w):
                    found = True
                    break
        if found:
            detected += 1

    return detected / n_trials


# ─── Experiment 1: mu-sweep ───────────────────────────────────────────────────

def experiment_mu_sweep(reference, p, w, L, depth, trials, delta,
                        mu_values, Q, n_avg, alpha):
    mu_c = mu_c_numerical(w, L, p, delta, depth)
    ws   = w_star(Q, n_avg, alpha)

    HR = '─' * 84
    print(f"\n{HR}")
    print(f"  EXPERIMENT 1: μ-sweep")
    print(f"  p={p:.1e}  w={w}  w*={ws}  L={L}  T={depth:,}  δ={delta}")
    print(f"  μ_c = {mu_c:.4f}   Note: q_approx is unreliable for low μ (see paper §2)")
    print(HR)
    print(f"  {'μ':>6}  {'q_exact':>9}  {'q_approx':>9}  "
          f"{'T*(theory)':>12}  {'Det theory':>11}  {'Det emp':>9}  {'|Δ|':>7}  {'Tier':>5}")
    print(HR)

    results = []
    for mu in mu_values:
        qe   = q_recurrence(w, L, mu)
        qa   = q_poisson(w, L, mu)
        Tth  = T_star(p, mu, w, L, delta)
        det_th = detection_prob_theory(p, mu, w, L, depth)
        det_em = simulate_detection(reference, mu, p, w, L, depth, trials)
        gap    = abs(det_em - det_th)
        tier   = '1' if mu == 0.0 else ('2' if mu <= mu_c else '3')
        Tstr   = f"{Tth:,}" if Tth != float('inf') else "∞"
        print(f"  {mu:>6.3f}  {qe:>9.5f}  {qa:>9.5f}  "
              f"{Tstr:>12}  {det_th:>11.4f}  {det_em:>9.4f}  {gap:>7.4f}  {tier:>5}")
        results.append(dict(mu=mu, q_exact=qe, q_approx=qa, T_star=Tth,
                            det_theory=det_th, det_empirical=det_em,
                            gap=gap, tier=tier))
    print(HR)
    return results, mu_c


# ─── Experiment 2: w-sweep ────────────────────────────────────────────────────

def experiment_w_sweep(reference, p, mu, L, depth, trials, delta,
                       w_values, Q, n_avg, alpha):
    ws    = w_star(Q, n_avg, alpha)
    T_opt = T_star(p, mu, ws, L, delta)

    HR = '─' * 84
    print(f"\n{HR}")
    print(f"  EXPERIMENT 2: w-sweep  (cost of over-specification, Theorem 3.4)")
    print(f"  p={p:.1e}  μ={mu}  L={L}  T={depth:,}  δ={delta}")
    print(f"  w*={ws}   T*(w*)={T_opt:,}")
    print(HR)
    print(f"  {'w':>4}  {'q_exact':>9}  {'T*(theory)':>12}  "
          f"{'Ratio T/T*(w*)':>15}  {'Det theory':>11}  {'Det emp':>9}")
    print(HR)

    results = []
    for w in w_values:
        qe    = q_recurrence(w, L, mu)
        Tth   = T_star(p, mu, w, L, delta)
        det_th = detection_prob_theory(p, mu, w, L, depth)
        det_em = simulate_detection(reference, mu, p, w, L, depth, trials)
        ratio  = (Tth / T_opt) if T_opt < float('inf') else float('inf')
        Tstr   = f"{Tth:,}" if Tth != float('inf') else "∞"
        rstr   = f"{ratio:.3f}" if ratio != float('inf') else "∞"
        mark   = "  ← w*" if w == ws else ""
        print(f"  {w:>4}  {qe:>9.5f}  {Tstr:>12}  "
              f"{rstr:>15}  {det_th:>11.4f}  {det_em:>9.4f}{mark}")
        results.append(dict(w=w, q_exact=qe, T_star=Tth, ratio=ratio,
                            det_theory=det_th, det_empirical=det_em))
    print(HR)
    return results


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="mNGS Pathogen Detection Simulator — validates Theorems 3.1, 3.3 and 3.4",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--genome',  default='synthetic:5000')
    parser.add_argument('--p',       type=float, default=1e-4)
    parser.add_argument('--w',       type=int,   default=24)
    parser.add_argument('--L',       type=int,   default=150)
    parser.add_argument('--delta',   type=float, default=0.01)
    parser.add_argument('--depth',   type=int,   default=100_000)
    parser.add_argument('--trials',  type=int,   default=300)
    # Divergence held fixed for Experiment 2. The paper reports the w-sweep at
    # mu = 0.10, near mu_c = 0.1349, where over-specifying w costs the most.
    parser.add_argument('--w_sweep_mu', type=float, default=0.10)
    parser.add_argument('--Q',       type=int,   default=5_000)
    parser.add_argument('--n_avg',   type=int,   default=15_000)
    parser.add_argument('--alpha',   type=float, default=1e-6)
    parser.add_argument('--out',     default='data/results.csv')
    parser.add_argument('--seed',    type=int,   default=42)
    args = parser.parse_args()

    random.seed(args.seed)

    if args.genome.startswith('synthetic:'):
        N = int(args.genome.split(':')[1])
        print(f"Generating synthetic reference genome: {N:,} bp ...")
        reference = random_sequence(N)
    else:
        print(f"Loading genome from {args.genome} ...")
        reference = load_fasta(args.genome)
    print(f"Reference length: {len(reference):,} bp")

    ws = w_star(args.Q, args.n_avg, args.alpha)
    print(f"w* (Theorem 3.4): {ws}")

    mu_vals = [0.00, 0.02, 0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20]
    mu_results, mu_c = experiment_mu_sweep(
        reference, args.p, args.w, args.L, args.depth, args.trials,
        args.delta, mu_vals, args.Q, args.n_avg, args.alpha
    )

    w_vals = [ws + k for k in range(0, 12, 2)]
    w_results = experiment_w_sweep(
        reference, args.p, args.w_sweep_mu, args.L, args.depth, args.trials,
        args.delta, w_vals, args.Q, args.n_avg, args.alpha
    )

    with open(args.out, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            'experiment', 'param_name', 'param_value',
            'q_exact', 'q_approx', 'T_star',
            'det_theory', 'det_empirical', 'gap', 'ratio', 'tier',
            'fixed_mu', 'trials'
        ])
        writer.writeheader()
        for r in mu_results:
            writer.writerow({'experiment': 'mu_sweep', 'param_name': 'mu',
                             'param_value': r['mu'], 'q_exact': r['q_exact'],
                             'q_approx': r['q_approx'], 'T_star': r['T_star'],
                             'det_theory': r['det_theory'],
                             'det_empirical': r['det_empirical'],
                             'gap': r['gap'], 'ratio': '', 'tier': r['tier'],
                             'fixed_mu': '', 'trials': args.trials})
        for r in w_results:
            writer.writerow({'experiment': 'w_sweep', 'param_name': 'w',
                             'param_value': r['w'], 'q_exact': r['q_exact'],
                             'q_approx': '', 'T_star': r['T_star'],
                             'det_theory': r['det_theory'],
                             'det_empirical': r['det_empirical'],
                             'gap': '', 'ratio': r['ratio'], 'tier': '',
                             'fixed_mu': args.w_sweep_mu,
                             'trials': args.trials})

    print(f"\nResults → {args.out}")
    print("Done.")


if __name__ == '__main__':
    main()
