#!/usr/bin/env python3
"""
Paper Figure Generator
=======================
Produces three publication-quality figures for:
  "Minimum Sequencing Depth Guarantees for Reference-Based Detection
   of Divergent Pathogens in Clinical Metagenomic Sequencing"

Figure 1 — Detection rate vs divergence (μ-sweep)
  Theoretical curve (exact recurrence) vs Poisson approximation
  vs empirical simulation. Marks μ_c and Tier boundary.

Figure 2 — Optimal w* and cost of over-specification
  Left: q_exact(w) vs w for several μ values.
  Right: T*(w)/T*(w*) ratio vs w, showing cost diverges as μ → μ_c.

Figure 3 — Parameter phase diagram (p, μ) → log T*
  Heatmap of log₁₀ T* over the full parameter space.
  Contour lines at practical depth limits.
  Bold μ_c boundary curve.

Usage:
    python plot_figures.py                    # uses results.csv if present
    python plot_figures.py --csv results.csv  # explicit CSV path
    python plot_figures.py --no-empirical     # theory-only (no CSV needed)
"""

import argparse
import csv
import math
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import numpy as np

# ── Reuse math from simulator (inline copies so this file is standalone) ─────

def q_recurrence(w, L, mu):
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
            total += p ** n
        P[n] = total
    return 1.0 - P[L]


def q_poisson(w, L, mu):
    lam = (L - w + 1) * ((1.0 - mu) ** w)
    return 1.0 - math.exp(-lam)


def T_star(p, mu, w, L, delta):
    q = q_recurrence(w, L, mu)
    if q <= 0.0:
        return float('inf')
    pq = p * q
    if pq <= 0.0:
        return float('inf')
    log_term = math.log(1.0 - pq)
    if log_term >= 0:
        return 1
    return math.ceil(math.log(delta) / log_term)


def detection_prob(p, mu, w, L, T, use_exact=True):
    q = q_recurrence(w, L, mu) if use_exact else q_poisson(w, L, mu)
    pq = p * q
    if pq <= 0.0:
        return 0.0
    return 1.0 - (1.0 - pq) ** T


def mu_c_numerical(w, L, p, delta, T_max, tol=1e-7):
    if T_max <= 0:
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
    return math.ceil(math.log(Q * n_avg / alpha) / math.log(4))


# ── Shared style ──────────────────────────────────────────────────────────────

COLORS = {
    'exact':     '#2166ac',   # blue
    'poisson':   '#d73027',   # red
    'empirical': '#1a9641',   # green
    'tier2':     '#d1e5f0',
    'tier3':     '#fddbc7',
    'muc_line':  '#636363',
    'wstar_line':'#762a83',
}

FONT = {'family': 'serif', 'size': 11}
plt.rcParams.update({
    'font.family':       'serif',
    'font.size':         11,
    'axes.labelsize':    12,
    'axes.titlesize':    13,
    'legend.fontsize':   10,
    'figure.dpi':        150,
    'axes.spines.top':   False,
    'axes.spines.right': False,
})


# ── Figure 1: detection rate vs μ ─────────────────────────────────────────────

def figure1(p, w, L, delta, depth, mu_c,
            emp_mu=None, emp_det=None,
            outpath='figure1_detection_vs_mu.pdf'):
    """
    Detection probability vs divergence μ.
    Shows exact theory, Poisson approximation, and empirical estimates.
    """
    mu_dense = np.linspace(0.0, 0.22, 300)

    det_exact   = [detection_prob(p, float(mu), w, L, depth, use_exact=True)
                   for mu in mu_dense]
    det_poisson = [detection_prob(p, float(mu), w, L, depth, use_exact=False)
                   for mu in mu_dense]

    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Tier shading
    ax.axvspan(0,    mu_c, alpha=0.18, color=COLORS['tier2'], label='Tier 2')
    ax.axvspan(mu_c, 0.22, alpha=0.18, color=COLORS['tier3'], label='Tier 3')

    # μ_c line
    ax.axvline(mu_c, color=COLORS['muc_line'], lw=1.4, ls='--', zorder=3)
    ax.text(mu_c + 0.003, 0.08,
            f'$\\mu_c = {mu_c:.3f}$',
            color=COLORS['muc_line'], fontsize=10)

    # δ target line
    ax.axhline(1 - delta, color='#888', lw=0.9, ls=':', zorder=2)
    ax.text(0.001, 1 - delta + 0.01,
            f'$1-\\delta = {1-delta:.2f}$',
            color='#888', fontsize=9)

    # Theory curves
    ax.plot(mu_dense, det_exact,
            color=COLORS['exact'], lw=2.2, zorder=5,
            label='Theorem 1 (exact recurrence)')
    ax.plot(mu_dense, det_poisson,
            color=COLORS['poisson'], lw=1.8, ls='--', zorder=4,
            label='Theorem 1 (Poisson approx.)')

    # Empirical points
    if emp_mu is not None and emp_det is not None:
        ax.scatter(emp_mu, emp_det,
                   color=COLORS['empirical'], s=45, zorder=6,
                   marker='o', label='Empirical (simulation)')

    ax.set_xlim(0, 0.22)
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel('Per-base divergence  $\\mu$')
    ax.set_ylabel('Detection probability')
    ax.set_title(
        f'Figure 1 — Detection rate vs divergence\n'
        f'$p={p:.0e}$,  $w={w}$,  $L={L}$,  $T={depth:,}$,  $\\delta={delta}$'
    )

    # Tier labels
    ax.text(mu_c / 2, 0.92, 'Tier 2\n(detectable)',
            ha='center', va='top', fontsize=9.5, color='#2166ac')
    ax.text(mu_c + (0.22 - mu_c) / 2, 0.92, 'Tier 3\n(effectively novel)',
            ha='center', va='top', fontsize=9.5, color='#d73027')

    ax.legend(loc='lower left', framealpha=0.9)
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved → {outpath}")


# ── Figure 2: optimal w* and over-specification cost ─────────────────────────

def figure2(p, L, delta, depth, Q, n_avg, alpha,
            mu_highlight=None,
            outpath='figure2_wstar_tradeoff.pdf'):
    """
    Left:  q_exact(w) vs w for several μ values.
    Right: T*(w)/T*(w*) ratio vs w — cost of over-specification.
    """
    if mu_highlight is None:
        mu_highlight = [0.02, 0.05, 0.10, 0.13]

    ws     = w_star(Q, n_avg, alpha)
    w_vals = list(range(ws - 2, ws + 18))

    mu_colors = plt.cm.plasma(np.linspace(0.15, 0.85, len(mu_highlight)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))

    for mu, col in zip(mu_highlight, mu_colors):
        qs     = [q_recurrence(int(w), L, mu) for w in w_vals]
        T_opt  = T_star(p, mu, ws, L, delta)
        ratios = []
        for w in w_vals:
            Tw = T_star(p, mu, int(w), L, delta)
            ratios.append(Tw / T_opt if T_opt < 1e15 else float('nan'))

        lbl = f'$\\mu = {mu}$'
        ax1.plot(w_vals, qs,    color=col, lw=2,   label=lbl)
        ax2.plot(w_vals, ratios, color=col, lw=2,   label=lbl)

    for ax in (ax1, ax2):
        ax.axvline(ws, color=COLORS['wstar_line'], lw=1.5, ls='--', zorder=5)
        ax.text(ws + 0.3, ax.get_ylim()[1] if ax == ax1 else
                max(1, ax.get_ylim()[1]) * 0.95,
                f'$w^*={ws}$',
                color=COLORS['wstar_line'], fontsize=10)

    ax1.set_xlabel('Window length  $w$')
    ax1.set_ylabel('Informative-read probability  $q(w,L,\\mu)$')
    ax1.set_title('Sensitivity vs window length')
    ax1.set_ylim(0, 1.05)
    ax1.legend(loc='lower left', framealpha=0.9)

    ax2.set_xlabel('Window length  $w$')
    ax2.set_ylabel('Required depth ratio  $T^*(w) / T^*(w^*)$')
    ax2.set_title('Cost of over-specification (Theorem 3)')
    ax2.axhline(1.0, color='#aaa', lw=0.9, ls=':')
    ax2.legend(loc='upper left', framealpha=0.9)

    fig.suptitle(
        f'Figure 2 — Optimal $w^* = {ws}$  '
        f'($Q={Q:,}$, $n_{{\\rm avg}}={n_avg:,}$, $\\alpha={alpha:.0e}$)',
        fontsize=12
    )
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved → {outpath}")


# ── Figure 3: phase diagram (p, μ) ────────────────────────────────────────────

def figure3(w, L, delta, T_limits=None,
            outpath='figure3_phase_diagram.pdf'):
    """
    Heatmap of log₁₀ T*(p, μ) over the biologically relevant parameter space.
    Contour lines show practical sequencing depth budgets.
    Bold curve shows μ_c(p) for each depth limit.
    """
    if T_limits is None:
        T_limits = [1e5, 1e6, 1e7]

    mu_vals = np.linspace(0.0, 0.22, 120)
    p_log   = np.linspace(-7, -2, 100)
    p_vals  = 10.0 ** p_log

    # Build T* grid
    Z = np.zeros((len(p_vals), len(mu_vals)))
    for i, p in enumerate(p_vals):
        for j, mu in enumerate(mu_vals):
            Tv = T_star(p, float(mu), w, L, delta)
            Z[i, j] = math.log10(Tv) if Tv < 1e15 else 15.0

    fig, ax = plt.subplots(figsize=(7.5, 5))

    img = ax.pcolormesh(mu_vals, p_log, Z,
                        cmap='YlOrRd', vmin=3, vmax=10, shading='auto')
    cbar = fig.colorbar(img, ax=ax, pad=0.02)
    cbar.set_label('$\\log_{10} T^*(p,\\mu,w,L,\\delta)$')
    cbar.set_ticks([3, 4, 5, 6, 7, 8, 9, 10])

    # Contour lines at practical depth limits
    limit_colors = ['#2166ac', '#4dac26', '#d01c8b']
    CS = ax.contour(mu_vals, p_log, Z,
                    levels=[math.log10(t) for t in T_limits],
                    colors=limit_colors, linewidths=1.8)
    labels = [f'$T={t:.0e}$' for t in T_limits]
    for i, (lbl, col) in enumerate(zip(labels, limit_colors)):
        CS.legend_elements()[0][i].set_label(lbl) if hasattr(CS, 'legend_elements') else None
    # Proxy artists for legend
    proxy = [plt.Line2D([0], [0], color=col, lw=1.8, label=lbl)
             for col, lbl in zip(limit_colors, labels)]

    # μ_c boundary for each depth limit
    for T_max, col in zip(T_limits, limit_colors):
        muc_vals = []
        for p in p_vals:
            muc_vals.append(mu_c_numerical(w, L, p, delta, T_max))
        ax.plot(muc_vals, p_log,
                color=col, lw=2.5, ls=':', zorder=6)

    ax.set_xlabel('Per-base divergence  $\\mu$')
    ax.set_ylabel('$\\log_{10}$(pathogen abundance  $p$)')
    ax.set_title(
        f'Figure 3 — Phase diagram: required sequencing depth\n'
        f'$w={w}$,  $L={L}$,  $\\delta={delta}$   '
        f'(dotted = $\\mu_c$ boundary for each depth limit)'
    )

    # Annotate Tier boundary qualitatively
    ax.text(0.21, -2.3, 'Tier 3\n(harder)', ha='right',
            fontsize=9, color='#d73027')
    ax.text(0.01, -6.7, 'Tier 2\n(accessible)', ha='left',
            fontsize=9, color='#2166ac')

    ax.legend(handles=proxy, loc='upper left', framealpha=0.9,
              title='Depth contours', title_fontsize=9)
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved → {outpath}")


# ── Load empirical data from CSV ──────────────────────────────────────────────

def load_empirical(csv_path):
    emp_mu, emp_det = [], []
    try:
        with open(csv_path) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                if row['experiment'] == 'mu_sweep':
                    emp_mu.append(float(row['param_value']))
                    emp_det.append(float(row['det_empirical']))
    except FileNotFoundError:
        print(f"  (CSV {csv_path} not found — skipping empirical overlay)")
    return (emp_mu, emp_det) if emp_mu else (None, None)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate three paper figures",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--csv',          default='results.csv')
    parser.add_argument('--no-empirical', action='store_true')
    parser.add_argument('--p',     type=float, default=1e-4)
    parser.add_argument('--w',     type=int,   default=24)
    parser.add_argument('--L',     type=int,   default=150)
    parser.add_argument('--delta', type=float, default=0.01)
    parser.add_argument('--depth', type=int,   default=100_000)
    parser.add_argument('--Q',     type=int,   default=5_000)
    parser.add_argument('--n_avg', type=int,   default=15_000)
    parser.add_argument('--alpha', type=float, default=1e-6)
    parser.add_argument('--outdir',        default='/mnt/user-data/outputs')
    args = parser.parse_args()

    import os
    os.makedirs(args.outdir, exist_ok=True)
    def out(name):
        return os.path.join(args.outdir, name)

    mu_c = mu_c_numerical(args.w, args.L, args.p, args.delta, args.depth)
    print(f"Parameters: p={args.p:.1e}  w={args.w}  L={args.L}  "
          f"δ={args.delta}  T={args.depth:,}")
    print(f"Computed μ_c = {mu_c:.4f}  |  "
          f"w* = {w_star(args.Q, args.n_avg, args.alpha)}")

    emp_mu, emp_det = (None, None) if args.no_empirical else \
        load_empirical(args.csv)

    print("\nGenerating figures ...")

    figure1(
        p=args.p, w=args.w, L=args.L, delta=args.delta,
        depth=args.depth, mu_c=mu_c,
        emp_mu=emp_mu, emp_det=emp_det,
        outpath=out('figure1_detection_vs_mu.pdf')
    )

    figure2(
        p=args.p, L=args.L, delta=args.delta, depth=args.depth,
        Q=args.Q, n_avg=args.n_avg, alpha=args.alpha,
        mu_highlight=[0.02, 0.05, 0.10, 0.13],
        outpath=out('figure2_wstar_tradeoff.pdf')
    )

    figure3(
        w=args.w, L=args.L, delta=args.delta,
        T_limits=[1e5, 1e6, 1e7],
        outpath=out('figure3_phase_diagram.pdf')
    )

    print("\nAll three figures saved.")


if __name__ == '__main__':
    main()
