# Minimum Sequencing Depth Guarantees for Reference-Based Detection of Divergent Pathogens in Clinical Metagenomic Sequencing

**Author:** Devansh Malhotra  
**Status:** Preprint on Figshare — [doi.org/10.6084/m9.figshare.33668701](https://doi.org/10.6084/m9.figshare.33668701)

---

## Abstract

Clinical metagenomic next-generation sequencing (mNGS) enables rapid identification of infectious agents by matching sequenced DNA reads from patient samples against a curated reference database. This approach is particularly valuable in time-critical settings where conventional culture-based diagnosis is too slow, yet it is most likely to fail in precisely the scenarios of greatest clinical concern: when a pathogen has diverged substantially from any entry in the reference database, as occurs with emerging variants, cross-species spillovers, and novel outbreak strains. In the absence of a principled mathematical framework, the minimum sequencing depth required to guarantee detection of a divergent pathogen remains determined empirically on a per-pathogen, per-assay basis, providing no general design criterion.

Here, a formal probabilistic model of reference-based $w$-mer classification is developed, and three theorems are established. A minimum sequencing depth bound $T^{\ast}(p, \mu, w, L, \delta)$ is derived as an explicit closed-form function of pathogen abundance $p$, per-base genomic divergence $\mu$, classifier window length $w$, read length $L$, and target false-negative rate $\delta$. The derivation rests on the exact probability $q(w, L, \mu)$ that a single read from a divergent pathogen contains at least one classifiable $w$-mer, computed via a linear recurrence on the run-length distribution of Bernoulli trials. The critical divergence $\mu_c$ (the maximum divergence at which detection can be guaranteed within a practical sequencing budget) is proved unique and computable by bisection on the exact recurrence, providing a quantitative boundary between classifiable and effectively novel pathogens. The optimal classifier window length $w^{\ast} = \lceil \log_4(Q\, n_{\mathrm{avg}} / \alpha) \rceil$ is shown to minimize required depth subject to a false-positive rate constraint, and the depth penalty for a longer window grows rapidly with divergence: a window 10 bases longer than $w^{\ast}$ requires $2.3\times$ more reads at $\mu = 0.10$ and $3.9\times$ more at $\mu_c$. A methodologically significant finding is that the Poisson approximation commonly used to estimate classifiable-read probability overestimates $q(w, L, \mu)$ by 72% at $\mu = 0.12$ alone, with the overestimate exceeding 100% as $\mu$ approaches $\mu_c$ from below, making the exact recurrence essential for conservative depth planning. The same approximation also overestimates the critical divergence itself (0.199 versus the exact 0.1349 at $T_{\max} = 10^5$). All results are validated by simulation against real SARS-CoV-2 (Wuhan) and Influenza A H3N2 HA reference genomes across biologically motivated divergence scenarios, with theory–empirical gaps of at most 0.062 across all experiments (0.021 on the real genomes); the residual gaps are themselves predicted by a shared-genome correction for reads drawn from a single mutated genome. The framework extends Ebinger et al. [1] by introducing the divergence dimension absent from their model and provides the first principled basis for sequencing depth selection in mNGS assay design when pathogen identity and divergence are uncertain.

> 📖 New to k-mer classification, divergence, or the math notation here?
> See [EXPLAINER.md](EXPLAINER.md) for a full plain-language walkthrough
> of every term, equation, and proof step.

---


## Repository Structure

```
├── code/
│   ├── simulator.py          # Component 1: Read simulator — reproduces paper Tables 1–2
│   ├── real_genomes.py       # Component 2: Real-genome validation — reproduces paper Table 3
│   ├── plot_figures.py       # Generates all three paper figures
│   ├── shared_genome_expectation.py  # Expected detection when all reads in a trial share one mutated genome (paper Sec. 5.4)
│   └── requirements.txt
├── proofs/
│   └── paper.tex             # Full LaTeX source with formal proofs
├── figures/
│   ├── fig1.png              # Detection rate vs divergence
│   ├── fig2.png              # Optimal w* and over-specification cost
│   └── fig3.png              # Phase diagram (p, mu) → T*
└── data/
    ├── sars_cov2_wuhan.fasta # Reference genome (MN908947, 29,903 bp)
    ├── h3n2_ha.fasta         # Influenza H3N2 HA (CY163680.1, 1,737 bp)
    ├── h1n1_ha.fasta         # Influenza H1N1 HA (CY121680.1, 1,752 bp)
    ├── results.csv           # Synthetic genome sweep results (regenerate: simulator.py)
    └── results_real_genomes.csv  # Real genome validation (regenerate: real_genomes.py)
```

---

## Quick Start

```bash
pip install -r code/requirements.txt

# Synthetic-genome sweeps -> data/results.csv (paper Tables 1-2)
python code/simulator.py

# Real-genome validation -> data/results_real_genomes.csv (paper Table 3)
python code/real_genomes.py

# All three figures, as both PDF and the PNGs embedded below -> figures/
python code/plot_figures.py

# Shared-genome correction quoted in paper Section 5.4
python code/shared_genome_expectation.py

# The same simulator on a different reference, with custom parameters
python code/simulator.py --genome data/sars_cov2_wuhan.fasta \
    --p 1e-4 --w 24 --L 150 --delta 0.01 \
    --depth 100000 --trials 200
```

---

## Research Question

> Given a clinical sample where a pathogen is present at relative abundance $p$ with per-base genomic divergence $\mu$ from the nearest reference in a database of size $Q$, what is the **minimum sequencing depth** $T$ that provably guarantees detection with false-negative rate $\leq \delta$ — as a closed-form function of $p$, $\mu$, classifier $w$-mer length $w$, and read length $L$?

Three sub-questions:

1. Can we derive a tight lower bound $T^{\ast}(p, \mu, w, L, \delta)$?
2. Where is the hard boundary $\mu_c$ beyond which no reference-based classifier can guarantee detection?
3. What is the optimal $w^{\ast}$ minimizing $T^{\ast}$ subject to a false-positive constraint?

---

## Model

### Setting

A clinical mNGS experiment sequences $T$ reads from a patient sample. A pathogen $\mathcal{P}$ is present at relative abundance $p \in (0,1)$, meaning fraction $p$ of reads originate from $\mathcal{P}$. The reference database $\mathcal{D}$ contains $Q$ sequences of average length $n_{\mathrm{avg}}$ bp; the nearest reference diverges from the true pathogen genome at per-base rate $\mu$ under an independent substitution model. $w$-mer classification follows the lowest-common-ancestor approach of Kraken2 [2].

### Definition 1.1 — Divergence Tier

Let $\mu \in [0,1]$ be the per-base divergence between pathogen and its nearest reference.

| Tier | Condition | Meaning |
|------|-----------|---------|
| **1 (Known)** | $\mu = 0$ | Exact reference in database |
| **2 (Divergent–detectable)** | $0 < \mu \leq \mu_c$ | Detection guaranteed within budget $T_{\max}$ |
| **3 (Effectively novel)** | $\mu > \mu_c$ | No reference-based classifier guarantees detection within $T_{\max}$ |

### Definition 1.2 — Informative Read

A read of length $L$ is **informative** if it contains at least one $w$-mer that exactly matches the corresponding reference window.

Under the independent-substitution model (each base matches reference with probability $1-\mu$), define:

$$q(w, L, \mu) = 1 - P_L(w, 1-\mu)$$

where $P_n = P_n(w, 1-\mu)$ is the probability that $n$ independent Bernoulli $(1-\mu)$ trials contain **no run of $w$ consecutive successes**. With the boundary condition

$$
P_n = 1 \qquad \text{for } 0 \leq n \leq w-1
$$

(no run of length $w$ fits in fewer than $w$ trials), $P_n$ satisfies the linear recurrence

$$
P_n = \sum_{j=0}^{w-1} (1-\mu)^j \cdot \mu \cdot P_{n-j-1}, \qquad n \geq w
$$

> **The boundary condition is not optional.** Replacing it with $P_0 = 1$ alone and truncating the sum at $\min(n-1, w-1)$ — a natural-looking alternative — gives $P_1 = \mu$ instead of $P_1 = 1$ and corrupts every subsequent value: at $\mu = 0.12$ it returns $q = 0.949$ instead of $0.580$, and at $\mu = 0.20$ it returns $0.824$ instead of $0.119$. Dropping the all-successes case makes every $P_n$ too small, so the truncated form overestimates $q$ at every $\mu \in (0,1)$ and therefore always *under*estimates $T^{\ast}$ — it fails dangerously, not safely. See paper Remark 2.5.

The **Poisson approximation** — which treats the $L-w+1$ windows as independent (cf. [3] for related k-mer frequency approximations) and is used here only for comparison:

$$
q(w, L, \mu) \approx 1 - \exp\bigl(-(L-w+1)(1-\mu)^{w}\bigr)
$$

> **Important:** The Poisson approximation overestimates $q$ by **72%** at $\mu = 0.12$ and by more than 100% as $\mu$ approaches $\mu_c$ (ratio $1.34$ at $\mu = 0.10$, $3.8$ at $\mu = 0.20$; negligible for $\mu \leq 0.05$). Matching windows occur in overlapping clusters of mean size roughly $1/\mu$, which the Poisson count ignores (see Figure 1). The exact recurrence is required.

### Definition 1.3 — Detection Event

Given $T$ reads, the false-negative event $\overline{\mathcal{E}}$ occurs if no informative pathogen read appears. The classifier **guarantees $\delta$-detection** if $\Pr(\overline{\mathcal{E}}) \leq \delta$.

---

## Main Results

### Theorem 1 — Minimum Sequencing Depth

*(Theorem 3.1 in the paper; the Corollary below is Corollary 3.2, Theorem 2 is Theorem 3.3, and Theorem 3 is Theorem 3.4.)*

**Statement.** The minimum sequencing depth guaranteeing $\delta$-detection is:

$$
T^{\ast}(p, \mu, w, L, \delta) = \left\lceil \frac{\ln(1/\delta)}{-\ln\bigl(1 - p \cdot q(w,L,\mu)\bigr)} \right\rceil \leq \left\lceil \frac{\ln(1/\delta)}{p \cdot q(w,L,\mu)} \right\rceil
$$

The right-hand expression is the small-$pq$ approximation: an **upper** bound on $T^{\ast}$, tight within a factor of 2 (up to rounding) whenever $pq \leq 0.5$.

**Proof.**
Each read is independently a pathogen read with probability $p$, and each pathogen read is independently informative with probability $q = q(w,L,\mu)$. Therefore each read is an informative pathogen read with probability $pq$, independently across reads.

The false-negative event requires all $T$ reads to be non-informative:

$$\Pr(\overline{\mathcal{E}}) = (1 - pq)^T$$

Setting $\Pr(\overline{\mathcal{E}}) \leq \delta$ and solving for $T$:

$$(1-pq)^T \leq \delta \iff T \ln(1-pq) \leq \ln\delta \iff T \geq \frac{\ln(1/\delta)}{-\ln(1-pq)}$$

Taking the ceiling gives the exact integer minimum. The inequality $-\ln(1-x) \geq x$ (equality only at $x=0$) gives $-\ln(1-pq) \geq pq$, hence $\ln(1/\delta)/(-\ln(1-pq)) \leq \ln(1/\delta)/(pq)$, which gives the right-hand bound. The factor-of-2 claim follows from $-\ln(1-x) \leq 2x$ for $x \leq 1/2$.

**Corollary (Monotonicity).** $T^{\ast}$ is strictly increasing in $\mu$ and strictly increasing in $w$, for all $\mu \in (0,1)$, $w \geq 1$, $L > w$.

*Proof.* $q(w,L,\mu) = 1 - P_L(w,1-\mu)$ is strictly decreasing in $\mu$ (larger divergence makes each base more likely to mismatch, reducing the probability of a run of $w$ consecutive matches) and strictly decreasing in $w$ (a longer run is harder to achieve for any $\mu > 0$). Since $T^{\ast}$ is strictly decreasing in $q > 0$, strict monotonicity in both $\mu$ and $w$ follows.

---

### Theorem 2 — Critical Divergence $\mu_c$ (paper Theorem 3.3)

**Statement.** Fix $w, L, p, \delta$ and an integer depth ceiling $T_{\max}$ with $p > 1-\delta^{1/T_{\max}}$ (so that a perfectly matching pathogen is detectable within $T_{\max}$). Let

$$\theta = \frac{1-\delta^{1/T_{\max}}}{p} \in (0,1),$$

which agrees with $\ln(1/\delta)/(pT_{\max})$ to relative error $\ln(1/\delta)/(2T_{\max})$. Define:

$$\mu_c(w, L, p, \delta, T_{\max}) = \sup\bigl\lbrace \mu \in [0,1] : T^{\ast}(p,\mu,w,L,\delta) \leq T_{\max} \bigr\rbrace$$

Then $\mu_c$ is the unique solution of $q(w, L, \mu_c) = \theta$, where $q$ is the **exact** informative-read probability, and $T^{\ast} \leq T_{\max}$ if and only if $\mu \leq \mu_c$. For $\mu > \mu_c$, no reference-based $w$-mer classifier achieves $\delta$-detection within $T_{\max}$ reads (Tier 3).

**Proof.**
Since $T_{\max}$ is an integer, $T^{\ast} \leq T_{\max}$ iff $\ln(1/\delta)/(-\ln(1-pq)) \leq T_{\max}$ iff $1-pq \leq \delta^{1/T_{\max}}$ iff $q \geq \theta$. The function $q(w,L,\cdot)$ is a polynomial in $\mu$ (by the recurrence), hence continuous on $[0,1]$, with $q(w,L,0) = 1 > \theta$ and $q(w,L,1) = 0 < \theta$. It is strictly decreasing on $(0,1)$ (by the Corollary), hence on $[0,1]$ by continuity. The intermediate value theorem and strict monotonicity give a unique $\mu_c$ with $q(w,L,\mu_c) = \theta$, and $q \geq \theta$ exactly when $\mu \leq \mu_c$.

**Computing $\mu_c$.** Solve $q(w,L,\mu_c) = \theta$ by bisection ($O(wL)$ per evaluation of $q$). Substituting the Poisson approximation gives a closed form,

$$\mu_c^{\mathrm{Pois}} = 1 - \left(\frac{-\ln(1-\theta)}{L-w+1}\right)^{1/w},$$

but it overestimates $\mu_c$ because the Poisson approximation overestimates $q$. For $p = 10^{-4}$, $w = 24$, $L = 150$, $\delta = 0.01$, $T_{\max} = 10^5$ the exact value is $\mu_c = 0.1349$, while the closed form gives $0.199$ (about 48% too high), which would misclassify $\mu = 0.15$ and $\mu = 0.18$ as Tier 2. All $\mu_c$ values reported here and in the paper are exact.

**Remark.** $\mu_c$ is not a fixed biological constant — it depends on all detection parameters. Increasing read length $L$ raises $\mu_c$ for all other fixed parameters: a formal justification for the practical preference for long-read platforms in novel-pathogen detection.

---

### Theorem 3 — Optimal Classifier Window Length $w^{\ast}$ (paper Theorem 3.4)

**Statement.** Let $\alpha \in (0,1)$ be a target false-positive rate. The false-positive probability for a single $w$-mer query against the database is:

$$
P_{\mathrm{FP}}(w) \approx \frac{Q \cdot n_{\mathrm{avg}}}{4^w}
$$

The optimal window length minimizing $T^{\ast}$ subject to $P_{\mathrm{FP}}(w) \leq \alpha$ is:

$$
w^{\ast} = \left\lceil \log_4\left(\frac{Q \cdot n_{\mathrm{avg}}}{\alpha}\right) \right\rceil
$$

For $\mu > 0$, using any $w > w^{\ast}$ strictly increases required sequencing depth, and the penalty grows rapidly with $\mu$.

**Proof.**
By the Corollary to Theorem 1, $T^{\ast}(w)$ is strictly increasing in $w$. Therefore the minimum over the feasible set $\mathcal{W} = \lbrace w \in \mathbb{Z}^+ : P_{\mathrm{FP}}(w) \leq \alpha \rbrace$ is attained at $w^{\dagger} = \min \mathcal{W}$.

The constraint $P_{\mathrm{FP}}(w) \leq \alpha$, with the false-positive model $P_{\mathrm{FP}}(w) \approx Q \cdot n_{\mathrm{avg}} / 4^w$ above, becomes:

$$
\frac{Q \cdot n_{\mathrm{avg}}}{4^w} \leq \alpha \iff 4^w \geq \frac{Q \cdot n_{\mathrm{avg}}}{\alpha} \iff w \geq \log_4\left(\frac{Q \cdot n_{\mathrm{avg}}}{\alpha}\right)
$$

The smallest integer satisfying this is $w^{\dagger} = \lceil \log_4(Q \cdot n_{\mathrm{avg}} / \alpha) \rceil = w^{\ast}$.

**Cost of over-specification.** For $w = w^{\ast} + \Delta$ with $\Delta \geq 1$, strict monotonicity in $w$ gives $T^{\ast}(w^{\ast}+\Delta) > T^{\ast}(w^{\ast})$ for every $\mu > 0$. When $pq \ll 1$ the ratio is approximately $q(w^{\ast}, L, \mu)/q(w^{\ast}+\Delta, L, \mu)$; for $w^{\ast} = 24$, $\Delta = 10$ it is $2.3$ at $\mu = 0.10$, $3.9$ at $\mu_c = 0.1349$, and $9.7$ at $\mu = 0.20$.

**Remark (Tradeoff sharpening under divergence).** At $\mu = 0$, $q(w,L,0) = 1$ for all $w \leq L$, so $T^{\ast}$ is independent of $w$ and any $w \geq w^{\ast}$ is acceptable. For $\mu > 0$, $T^{\ast}$ grows strictly in $w$ beyond $w^{\ast}$ and the numerical penalty grows rapidly with $\mu$ — when correct window choice matters most.

---

## Derivation of the Linear Recurrence (Eq. 1)

For $n < w$ no run of length $w$ fits in the sequence at all, which gives the boundary condition $P_n = 1$. For $n \geq w$, the recurrence is derived by conditioning on the position of the **first mismatch** in the Bernoulli sequence. If the first mismatch occurs at position $j+1$ for $j \in \lbrace 0, \ldots, w-1 \rbrace$:

- Positions $1, \ldots, j$ are all successes (matches), probability $(1-\mu)^j$
- Position $j+1$ is a mismatch, probability $\mu$
- No run of $w$ consecutive matches has been created (since $j < w$)
- The remaining $n-j-1$ positions must have no such run: $P_{n-j-1}$

If instead no mismatch occurs in the first $w$ trials, a run of length $w$ is already complete — an outcome excluded from $P_n$. These $w+1$ cases are exhaustive and mutually exclusive, so summing over $j$ yields the recurrence. As a check, the coefficients sum to $\sum_{j=0}^{w-1}(1-\mu)^j\mu = 1-(1-\mu)^w$, the probability that the first $w$ trials are not all successes. The recurrence has **order $w$** and evaluates in $O(wL)$ time.

This is equivalent to the classical "run probability" problem in combinatorics — the probability of obtaining at least one run of $w$ consecutive successes in $L$ Bernoulli trials — appearing here for the first time in a bioinformatics detection-guarantee context. The physical genome coverage analogue (read depth vs. assembly completeness) was first treated by Lander and Waterman [4]; Theorem 1 here is the classifier-detection counterpart of that foundational result.

---

## Numerical Calibration

With $Q = 5{,}000$, $n_{\mathrm{avg}} = 15{,}000$ bp, $\alpha = 10^{-6}$:

$$
w^{\ast} = \left\lceil \log_4\left(\frac{5000 \times 15000}{10^{-6}}\right) \right\rceil = \left\lceil \log_4(7.5 \times 10^{13}) \right\rceil = \lceil 23.05 \rceil = 24
$$

At $p = 10^{-4}$, $w = 24$, $L = 150$, $\delta = 0.01$, $T_{\max} = 100{,}000$: $\mu_c = 0.1349$

| Parameter set | $T^{\ast}$ (theory) | Det. rate (emp.) |
|---|---|---|
| $\mu = 0.00$ (Tier 1) | 46,050 | 1.0000 |
| $\mu = 0.05$ (Tier 2) | 46,705 | 1.0000 |
| $\mu = 0.12$ (Tier 2, near boundary) | 79,409 | 0.9967 |
| $\mu = 0.15$ (Tier 3) | 130,987 | 0.9700 |
| $\mu = 0.20$ (Tier 3) | 387,793 | 0.6333 |

Empirical rates: synthetic 5,000 bp genome, 300 trials per row, seed 42 (paper Table 1).

---

## Real Genome Validation

Reference genomes were retrieved from Nextstrain [5] and NCBI. Divergence values for SARS-CoV-2 variants are based on published SNP counts [6,7]; influenza inter-season drift estimates follow Bedford et al. [8].

Each row mutates the **reference** genome named below at the stated $\mu$ and queries the result against that same reference. The two cross-family rows therefore model a divergent pathogen rather than aligning two downloaded genomes: the accession in parentheses is the organism the divergence value is calibrated to, not a second reference. This is why the cross-subtype row reports $N = 1{,}737$ (H3N2 HA) and not 1,752 (H1N1 HA).

| Scenario | Reference | $N$ (bp) | $\mu$ | $q_{\mathrm{exact}}$ | Det. theory | Det. empirical | SE | $\|\Delta\|$ |
|---|---|---|---|---|---|---|---|---|
| SARS-CoV-2 Wuhan → Wuhan | MN908947 [6] | 29,903 | 0.0000 | 1.00000 | 1.0000 | 1.0000 | 0.0000 | 0.0000 |
| SARS-CoV-2 Wuhan → Omicron BA.1 | MN908947 (vs OX049951 [7]) | 29,903 | 0.0017 | 1.00000 | 1.0000 | 1.0000 | 0.0000 | 0.0000 |
| SARS-CoV-2 Wuhan → SARS-CoV-1 | MN908947 (vs AY278741) | 29,903 | 0.2000 | 0.11875 | 0.6950 | 0.6840 | 0.0104 | 0.0110 |
| H3N2 HA: same season | CY163680 [5] | 1,737 | 0.0000 | 1.00000 | 1.0000 | 1.0000 | 0.0000 | 0.0000 |
| H3N2 HA: 1-year drift | CY163680 [5] | 1,737 | 0.0100 | 1.00000 | 1.0000 | 0.9995 | 0.0005 | 0.0005 |
| H3N2 HA: 5-year drift | CY163680 [5] | 1,737 | 0.0500 | 0.98598 | 0.9999 | 1.0000 | 0.0000 | 0.0001 |
| H3N2 vs H1N1 HA (cross-subtype) | CY163680 (vs CY121680 [5]) | 1,737 | 0.3000 | 0.00742 | 0.0715 | 0.0505 | 0.0049 | 0.0210 |

2,000 trials per row, seed 42. Reproduce with `python code/real_genomes.py`. SE is the binomial standard error of the Monte Carlo estimate.

Maximum theory–empirical gap on the real genomes: **0.0210**. The largest gap anywhere in the paper is 0.062 (synthetic genome, $\mu = 0.20$, Table 1).

Both Tier 3 rows fall *below* theory, and the size of each shortfall is predicted quantitatively by the shared-genome correction (paper Section 5.4), which accounts for all reads in a trial coming from one mutated genome. Evaluating that correction independently of the run above gives expected detection rates of **0.6895** and **0.0517**, against the observed $0.6840 \pm 0.0104$ and $0.0505 \pm 0.0049$ — within one standard error in both cases. Run `python code/shared_genome_expectation.py` to reproduce. The gap is largest for the cross-subtype row because it combines the shortest genome with the highest divergence; the bias has a known sign, so $T^{\ast}$ stays conservative.

---

## Figures

### Figure 1 — Detection Rate vs Divergence

![Figure 1](figures/fig1.png)

Exact-recurrence theory (blue) tracks the empirical rate closely for $\mu \leq 0.15$; at $\mu = 0.18$–$0.20$ empirical detection falls modestly below theory (largest gap 0.062 in the paper's Table 1), as expected when all reads in a trial come from one mutated genome (paper Section 5.4). The Poisson approximation (red dashed) overestimates detection once $\mu$ is moderate: a 72% overestimate of $q$ at $\mu = 0.12$, growing to $3.8\times$ at $\mu = 0.20$. $\mu_c = 0.1349$ (vertical dashed line) marks the Tier 2 / Tier 3 boundary at $T_{\max} = 10^5$.

### Figure 2 — Optimal $w^{\ast}$ and Cost of Over-Specification

![Figure 2](figures/fig2.png)

Left: $q(w,L,\mu)$ decreases strictly with $w$ (Corollary). Right: $T^{\ast}(w)/T^{\ast}(w^{\ast})$ rises from 1.0 at $w = w^{\ast} = 24$, with steeper cost at higher $\mu$ — confirming Theorem 3 (paper Theorem 3.4).

### Figure 3 — Parameter Phase Diagram

![Figure 3](figures/fig3.png)

Heatmap of $\log_{10} T^{\ast}(p, \mu)$ over the biologically relevant parameter space. Contour lines at $10^5$, $10^6$, $10^7$ reads correspond to practical sequencing platform limits. Dotted curves show $\mu_c(p)$ for each depth budget.

---

## Key Novel Contributions

1. $T^{\ast}(p, \mu, w, L, \delta)$ — first closed-form depth bound incorporating divergence
2. $\mu_c$ — first quantitative definition of the Tier 2 / Tier 3 boundary
3. $w^{\ast} = 24$ proved, not just empirical — first formal derivation of the optimal window
4. Poisson approximation failure — overestimates $q$ by 72% at $\mu = 0.12$ and by more than 100% near $\mu_c$, and overestimates $\mu_c$ itself ($0.199$ vs. exact $0.1349$); exact recurrence is necessary

---

## Relationship to Prior Work

| Prior work | What it covers | Gap addressed by this paper |
|---|---|---|
| Ebinger et al. (2021) [1] | $T \geq f(p, N, \delta)$ via Bernoulli | No divergence $\mu$; treats all reads as equally informative |
| Hera, Pierce-Ward, Koslicki (2023) [3] | k-mer statistics under mutation | For *estimating* divergence, not guaranteeing detection |
| Lander & Waterman (1988) [4] | Physical genome coverage for assembly | Not about classifier-based detection |
| Wood, Lu & Langmead (Kraken2, 2019) [2] | Empirical sensitivity per organism | No analytical false-negative bound |
| **This work** | $T^{\ast}(p, \mu, w, L, \delta)$ + $\mu_c$ + $w^{\ast}$ | Fills all three gaps simultaneously |

---

## Limitations

- **Independent substitution model:** Real divergence includes insertions, deletions, and recombination.
- **Read independence:** Theorem 1 treats pathogen reads as independent, but all pathogen reads in a sample derive from one mutated genome, so expected detection falls slightly below theory (Jensen's inequality), most noticeably for short, highly diverged genomes. The size of the shortfall is predicted to within one standard error for both Tier 3 rows of the real-genome table above. The bias has a known sign, so $T^{\ast}$ remains conservative. See paper Section 5.4 and `code/shared_genome_expectation.py`.
- **Uniform read sampling:** GC bias and amplification artifacts cause non-uniform coverage.
- **Sequencing errors:** Can be incorporated by replacing $\mu$ with $\mu_{\mathrm{eff}} = \mu + e - \mu e$ where $e$ is per-base error rate.
- **Tier 3 scope:** Reference-free detection of truly novel pathogens ($\mu > \mu_c$) is explicitly out of scope.

---

## Citation

```
Malhotra, D. (2026). Minimum Sequencing Depth Guarantees for Reference-Based
Detection of Divergent Pathogens in Clinical Metagenomic Sequencing.
Figshare preprint. https://doi.org/10.6084/m9.figshare.33668701
```

---

## References

[1] Ebinger, A., Fischer, S., & Höper, D. (2021). A theoretical and generalized approach for the assessment of the sample-specific limit of detection for clinical metagenomics. *Computational and Structural Biotechnology Journal*, 19, 732–742. https://doi.org/10.1016/j.csbj.2020.12.040

[2] Wood, D.E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20, 257. https://doi.org/10.1186/s13059-019-1891-0

[3] Hera, M.R., Pierce-Ward, N.T., & Koslicki, D. (2023). Deriving confidence intervals for mutation rates across a wide range of evolutionary distances using FracMinHash. *Genome Research*, 33(7), 1061–1068. https://doi.org/10.1101/gr.277651.123 (see also: Blanca, A., Harris, R.S., Koslicki, D., & Medvedev, P. (2022). The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches. *Journal of Computational Biology*, 29(2), 155–168. https://doi.org/10.1089/cmb.2021.0431)

[4] Lander, E.S., & Waterman, M.S. (1988). Genomic mapping by fingerprinting random clones: a mathematical analysis. *Genomics*, 2(3), 231–239. https://doi.org/10.1016/0888-7543(88)90007-9

[5] Nextstrain seasonal influenza repository. (2024). Reference sequences for H3N2 HA (A/Wisconsin/67/2005, CY163680.1) and H1N1pdm HA (A/California/07/2009, CY121680.1). https://github.com/nextstrain/seasonal-flu

[6] Wu, F., Zhao, S., Yu, B., et al. (2020). A new coronavirus associated with human respiratory disease in China. *Nature*, 579, 265–269. https://doi.org/10.1038/s41586-020-2008-3 (GenBank: MN908947)

[7] Viana, R., Moyo, S., Amoako, D.G., et al. (2022). Rapid epidemic expansion of the SARS-CoV-2 Omicron variant in southern Africa. *Nature*, 603, 679–686. https://doi.org/10.1038/s41586-022-04411-y (GenBank: OX049951 / BA.1 lineage divergence ~50 SNPs from MN908947)

[8] Bedford, T., Cobey, S., Beerli, P., & Pascual, M. (2010). Global migration dynamics underlie evolution and persistence of human influenza A (H3N2). *PLoS Pathogens*, 6(5), e1000918. https://doi.org/10.1371/journal.ppat.1000918

---

## License

GNU Affero General Public License v3.0 (AGPL-3.0). See [LICENSE](LICENSE) for full terms.

Reference genomes sourced from Nextstrain (nextstrain.org) under their respective open-data licenses.
