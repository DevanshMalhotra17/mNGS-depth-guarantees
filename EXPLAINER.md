# The Full Explainer: Everything You Need to Understand This Research

**Who this is for:** Anyone curious about this project: whether you're a biology student, a math student, a software developer, or just someone who wants to understand what's actually going on. No prior knowledge of bioinformatics, statistics, or advanced math is assumed. Everything is built from scratch.

**What this covers:** Every technical term, every equation, and every proof step explained as if you've never seen it before.

---

## Table of Contents

1. [The Big Picture: What Problem Are We Solving?](#1-the-big-picture)
2. [Background: How Clinical DNA Sequencing Works](#2-background-how-clinical-dna-sequencing-works)
3. [The Core Problem: Divergent Pathogens](#3-the-core-problem-divergent-pathogens)
4. [Setting Up the Math: Variables and What They Mean](#4-setting-up-the-math)
5. [Definition 1.1: The Three Tiers of Divergence](#5-definition-11-the-three-tiers)
6. [Definition 1.2: What Makes a Read "Informative"?](#6-definition-12-informative-reads)
7. [The Linear Recurrence: Exact Math Behind q](#7-the-linear-recurrence)
8. [Theorem 1: How Deep Do We Need to Sequence?](#8-theorem-1-minimum-depth) — paper Theorem 3.1
9. [Theorem 2: The Hard Boundary — Critical Divergence](#9-theorem-2-critical-divergence) — paper Theorem 3.3
10. [Theorem 3: Choosing the Right Window Size](#10-theorem-3-optimal-window-size) — paper Theorem 3.4
11. [Walk-Through: Real Numbers](#11-real-number-walk-through)
12. [Why the Poisson Approximation Fails](#12-why-the-poisson-approximation-fails)
13. [Glossary of Every Term Used](#13-glossary)

---

## 1. The Big Picture

Imagine a patient arrives in a hospital emergency room with a severe infection. Doctors need to know: **what pathogen is causing this?** Is it a bacteria? A virus? Which one? This matters enormously because the wrong antibiotic or the wrong treatment can kill the patient.

The traditional approach is to take a sample (blood, spinal fluid, etc.), grow it in a lab culture, and wait 24-72 hours for results. For critically ill patients, this is often too slow.

**Metagenomic next-generation sequencing (mNGS)** is a modern alternative: instead of growing anything, you sequence all the DNA in the sample directly (pathogen and human DNA together) and use a computer to figure out which pathogen is present. On long-read platforms this has been done in about 6 hours, start to finish.

But there's a catch: the computer works by comparing the sequenced DNA to a **reference database** of known pathogens. If the pathogen has **mutated significantly** from anything in the database or if it's a truly novel organism, the computer might miss it entirely.

**This paper answers the question:** Given that a pathogen has diverged (mutated) by some amount from the nearest known reference, how much sequencing do you need to do to be mathematically guaranteed to catch it?

This is not answered anywhere in the existing scientific literature. We provide the first closed-form (exact formula) answer.

---

## 2. Background: How Clinical DNA Sequencing Works

### What is DNA?

DNA is a molecule that carries genetic information. It's made of four chemical "letters" called **bases**: **A** (adenine), **T** (thymine), **G** (guanine), and **C** (cytosine). The entire genetic code of any organism is just a very long string of these four letters. A human genome is about 3 billion letters long. A coronavirus genome is about 30,000 letters long.

### What is sequencing?

**Sequencing** means reading the order of these letters. Modern sequencing machines (Illumina, Oxford Nanopore, etc.) can read millions of short snippets of DNA simultaneously.

### What is next-generation sequencing (NGS)?

NGS is a family of modern sequencing technologies that can read millions or billions of DNA fragments at once, very cheaply and quickly. The previous generation could only read one fragment at a time.

### What is metagenomic sequencing (mNGS)?

**Metagenomic** means you're sequencing ALL the DNA in a sample, without knowing in advance what's there. In a patient blood sample, most of the DNA will be human, but a small fraction might belong to a pathogen (a disease-causing organism). mNGS sequences everything and then sorts it out computationally.

### What is a read?

When a sequencing machine runs, it produces millions of short DNA snippets called **reads**. Each read is a fragment of the original DNA, typically 75–300 base pairs (letters) long. The machine essentially picks random positions in the DNA pool and reads short stretches from there.

### How does a computer identify a pathogen from reads?

It uses a **reference database**: a collection of known pathogen genome sequences. The computer compares each read to the database and asks: does any part of this read match a known pathogen?

The specific method used is called **k-mer classification**. A k-mer is just a substring of length k. For example, in the sequence `ATGCGA`, the 3-mers (k=3) are: `ATG`, `TGC`, `GCG`, `CGA`. A classifier takes each read, extracts all its k-mers, and checks how many of them appear in the reference database.

---

## 3. The Core Problem: Divergent Pathogens

### What is divergence?

Pathogens mutate over time. When a virus replicates, it occasionally makes copying errors. These errors called **mutations** change individual letters in the genome. Over many generations, enough mutations accumulate that the genome looks noticeably different from its ancestors.

**Divergence** ($\mu$ in our paper) is the fraction of letters that have changed. A divergence of $\mu = 0.05$ means 5% of the genome's letters are different from the reference. For a 30,000-letter virus, that's 1,500 different letters.

### Why does divergence cause detection failure?

K-mer classifiers rely on **exact matches**. A k-mer in a patient's read must match a k-mer in the reference database exactly — letter for letter. If even one letter in a k-mer is mutated, the match fails.

So if a pathogen has diverged significantly, its k-mers look different from the reference, and the classifier may not recognize any of them. The patient's infection goes undetected.

### Why is this a biosecurity problem?

In a pandemic or bioterrorism scenario, a novel or rapidly mutating pathogen might differ substantially from anything in the database. Without a mathematical guarantee, we don't know how much sequencing is needed before we can confidently say "we would have caught it if it were there."

---

## 4. Setting Up the Math

Before we write any equations, let's define every variable clearly.

### Variables and their meaning

| Symbol | Name | Plain English meaning |
|--------|------|----------------------|
| $p$ | Pathogen abundance | What fraction of all the DNA in the sample comes from the pathogen. If there are 1,000,000 total DNA molecules and 100 are from the pathogen, $p = 0.0001 = 10^{-4}$. |
| $\mu$ | Per-base divergence | For each letter in the pathogen genome, the probability that it differs from the nearest reference. $\mu = 0.05$ means 5% of letters have mutated. |
| $w$ | Window length (k-mer size) | The length of the substrings we're matching. We look for exact matches of length $w$. Also written as $k$ in other papers. |
| $L$ | Read length | How long each DNA snippet is that comes off the sequencing machine. Typically 75–300 base pairs. |
| $T$ | Sequencing depth | How many total reads the machine produces. More depth = more chances to see a rare pathogen. |
| $\delta$ | False-negative rate target | The maximum acceptable probability that we miss a pathogen that's actually there. $\delta = 0.01$ means we accept at most a 1% chance of missing it. |
| $Q$ | Database size | How many reference sequences are in our database. |
| $n_{\mathrm{avg}}$ | Average reference length | The average length (in base pairs) of sequences in the database. |
| $\alpha$ | False-positive rate target | The maximum acceptable probability that an innocent (non-pathogen) k-mer accidentally matches the database. |
| $N$ | Genome length | How long the pathogen's genome is in base pairs. |

### The core question as an equation

We want to find the minimum $T$ (number of reads) such that:

$$\text{Probability of missing the pathogen} \leq \delta$$

Our answer is a formula: $T^{\ast}(p, \mu, w, L, \delta)$ — the minimum depth needed, as a function of all the relevant parameters.

---

## 5. Definition 1.1: The Three Tiers

We divide all possible pathogens into three groups based on how diverged they are.

### Tier 1: Known ($\mu = 0$)

The pathogen matches the reference exactly. Every k-mer in every read will match the database. Detection is essentially guaranteed with any reasonable sequencing depth.

**Example:** SARS-CoV-2 Wuhan sequenced against the Wuhan reference genome.

### Tier 2: Divergent but Detectable ($0 < \mu \leq \mu_c$)

The pathogen has mutated, but not so much that detection becomes impossible. With enough sequencing, we can still guarantee we'd catch it. The exact boundary $\mu_c$ (called the "critical divergence") is calculated by Theorem 2.

**Example:** SARS-CoV-2 Omicron BA.1 has about 50 mutations compared to the Wuhan reference in a 29,903-letter genome. That's $\mu \approx 0.0017$ — well within Tier 2.

### Tier 3: Effectively Novel ($\mu > \mu_c$)

The pathogen has diverged so much that no matter how much you sequence, a reference-based k-mer classifier simply won't catch it. The k-mers look too foreign. This organism needs **reference-free methods** — algorithms that don't rely on having a matching example in the database.

**Example:** SARS-CoV-2 compared to SARS-CoV-1 (the 2003 SARS virus). These are related but differ by about 20% of their genome ($\mu \approx 0.20$), placing SARS-CoV-1 in Tier 3 when using a SARS-CoV-2 reference database.

---

## 6. Definition 1.2: Informative Reads

### What is an "informative" read?

A read is **informative** if the sequencing machine happened to capture a piece of the pathogen genome, AND that piece contains at least one substring of length $w$ that exactly matches the reference database.

Think of it like a puzzle: the read is a puzzle piece, and the reference database is the picture on the box. A read is informative if at least one small $w$-length section of it matches the picture perfectly.

### Why "at least one" match?

Because a read of length $L$ contains $L - w + 1$ overlapping windows of length $w$. You only need one of them to match. If the pathogen has diverged at rate $\mu$, some windows will be fully mutated (no match), but others might have all $w$ letters still matching the reference by chance.

### The formula for q(w, L, μ)

$q(w, L, \mu)$ is the probability that a single read from the divergent pathogen is informative.

We model each letter of the read as independently matching the reference with probability $1 - \mu$ (the probability of NOT being mutated). Think of each letter as a coin flip: heads (probability $1-\mu$) means it matches; tails (probability $\mu$) means it's mutated.

A window of length $w$ is an exact match only if **all $w$ letters** in that window are heads (matches). The probability of this for any specific window is:

$$(1-\mu)^w$$

For example: with $\mu = 0.05$ and $w = 24$:
$$(1 - 0.05)^{24} = 0.95^{24} \approx 0.292$$

So about 29% of windows would be fully matching.

### Why is the exact formula more complex?

Because windows in a read **overlap** — they share most of their letters. Window 1 covers positions 1–24, window 2 covers positions 2–25, and so on. These windows are not independent of each other; they share 23 letters. This correlation makes the calculation more involved than a simple probability calculation.

The **exact** probability that at least one window matches is:

$$q(w, L, \mu) = 1 - P_L(w,\, 1-\mu)$$

where $P_L(w, 1-\mu)$ is the probability that **no** window of length $w$ is fully matching in a sequence of $L$ letters. This is calculated via the linear recurrence in Section 7.

---

## 7. The Linear Recurrence

### What is a Bernoulli trial?

A **Bernoulli trial** is a single experiment with two outcomes: success (probability $p$) or failure (probability $1-p$). Flipping a coin is a Bernoulli trial (heads = success, tails = failure).

In our context: each letter in a read is a Bernoulli trial. Success (probability $1-\mu$) = the letter matches the reference. Failure (probability $\mu$) = the letter is mutated.

### What is a "run"?

A **run** of successes is a consecutive sequence of successful Bernoulli trials. In a coin flip sequence like H-H-T-H-H-H-T, the longest run of heads is 3.

In our context: a run of $w$ successes means $w$ consecutive letters in a read all match the reference — i.e., a complete window match.

### The question we're answering

What is the probability that in $L$ Bernoulli trials (letters), there is **no** run of $w$ consecutive successes (complete window match)?

Call this probability $P_n$ for $n$ trials (we'll eventually want $P_L$).

### Building the recurrence step by step

We use a clever conditioning trick. Instead of trying to compute the probability directly, we ask: **where does the first failure (mismatch) occur?**

Case 1: The first failure is at position 1 (the very first letter is a mismatch).
- This happens with probability $\mu$
- After position 1, we still have $n-1$ letters left, and they must contain no run of $w$ successes: probability $P_{n-1}$
- Contribution: $\mu \cdot P_{n-1}$

Case 2: The first failure is at position 2 (first letter matches, second doesn't).
- This happens with probability $(1-\mu) \cdot \mu$
- After position 2, we have $n-2$ letters left: probability $P_{n-2}$
- Contribution: $(1-\mu) \cdot \mu \cdot P_{n-2}$

Case 3: The first failure is at position 3 (first two match, third doesn't).
- Probability $(1-\mu)^2 \cdot \mu$, remaining: $P_{n-3}$
- Contribution: $(1-\mu)^2 \cdot \mu \cdot P_{n-3}$

...continuing this pattern until...

Case $w$: The first failure is at position $w$ (first $w-1$ letters match, the $w$-th doesn't).
- Probability $(1-\mu)^{w-1} \cdot \mu$, remaining: $P_{n-w}$
- Contribution: $(1-\mu)^{w-1} \cdot \mu \cdot P_{n-w}$

**But why stop at case $w$?** Because if the first failure is at position $w+1$ or later, that means positions 1 through $w$ are ALL successes — that's already a run of $w$ successes, which is exactly what we're trying to avoid. So those cases contribute 0 to $P_n$.

Adding up all the cases, for $n \geq w$:

$$P_n = \sum_{j=0}^{w-1} (1-\mu)^j \cdot \mu \cdot P_{n-j-1}, \qquad n \geq w$$

**Boundary condition:** $P_n = 1$ for every $n$ from $0$ to $w-1$.

$$P_n = 1 \qquad \text{for } 0 \leq n \leq w-1$$

**Why all of $0$ through $w-1$, and not just $P_0 = 1$?** Because a run of $w$ consecutive successes simply does not fit in fewer than $w$ trials, so for any such $n$ the probability of having no run is exactly 1 — not just for the empty sequence. This is the single easiest thing to get wrong in the whole derivation, and getting it wrong is not a rounding error. If you set only $P_0 = 1$ and let the sum run from $j=0$ to $\min(n-1, w-1)$ for all $n \geq 1$, you get $P_1 = \mu$ instead of $P_1 = 1$, and every value after that inherits the mistake: at $\mu = 0.12$ you would compute $q = 0.949$ instead of the correct $0.580$, and at $\mu = 0.20$ you would compute $0.824$ instead of $0.119$. Dropping the all-successes case makes every $P_n$ too small, so the truncated form overestimates $q$ at *every* $\mu$ in $(0,1)$ — and since $T^{\ast}$ falls as $q$ rises, it always tells you that you need **fewer** reads than you really do. The error fails dangerously, not safely. (Paper Remark 2.5 works this through.)

A useful sanity check on the recurrence itself: the coefficients sum to $\sum_{j=0}^{w-1}(1-\mu)^j \mu = 1 - (1-\mu)^w$, which is exactly the probability that the first $w$ trials are *not* all successes. That has to be true, because the $j$ cases above are precisely the ways the first $w$ trials can fail to be all-successes.

### Worked example with small numbers

Let's compute $P_n$ for $w = 2$ (we want no run of 2 consecutive matches), $\mu = 0.2$, so $1-\mu = 0.8$.

$P_0 = 1$ and $P_1 = 1$: both come straight from the boundary condition, since $n = 0$ and $n = 1$ are both below $w = 2$ and you cannot fit a run of 2 into fewer than 2 trials. Notice what the broken version would give here: summing $j=0$ alone yields $P_1 = 0.2$, which is already visibly wrong — one coin flip cannot produce two consecutive heads, so the answer has to be 1.

$P_2$: now $n \geq w$, so the recurrence applies. Two letters, which can have a run of 2 only if both match.
- $j=0$: $0.2 \cdot P_1 = 0.2 \cdot 1.0 = 0.20$
- $j=1$: $0.8 \cdot 0.2 \cdot P_0 = 0.8 \cdot 0.2 \cdot 1 = 0.16$
- $P_2 = 0.20 + 0.16 = 0.36$

**Check:** $P(\text{both match}) = 0.8^2 = 0.64$. So $P(\text{no run of 2}) = 1 - 0.64 = 0.36$. ✓

$P_3$: Three letters.
- $j=0$: $0.2 \cdot P_2 = 0.2 \cdot 0.36 = 0.072$
- $j=1$: $0.8 \cdot 0.2 \cdot P_1 = 0.16 \cdot 1.0 = 0.160$
- $P_3 = 0.072 + 0.160 = 0.232$

**Check by enumeration:** All length-3 sequences with no "SS" (two consecutive successes):
- FFF: $0.2^3 = 0.008$
- FFS: $0.2^2 \times 0.8 = 0.032$
- FSF: $0.2 \times 0.8 \times 0.2 = 0.032$
- SFF: $0.8 \times 0.2^2 = 0.032$
- SFS: $0.8 \times 0.2 \times 0.8 = 0.128$
- Total: $0.008 + 0.032 + 0.032 + 0.032 + 0.128 = 0.232$ ✓

The recurrence is verified correct.

### Computing q

Once we compute $P_L$ using the recurrence, the informative-read probability is simply:

$$q(w, L, \mu) = 1 - P_L$$

---

## 8. Theorem 1: Minimum Depth

*This is **Theorem 3.1** in the paper (`proofs/paper.tex`), Eq. 6. The "Theorem 1/2/3" numbering used in this explainer is local to this document; the paper numbers theorems by section.*

### The setup

We have $T$ reads coming off the sequencing machine. Each read is, independently:
- A pathogen read with probability $p$ (if the pathogen makes up fraction $p$ of the DNA pool)
- And if it's a pathogen read, informative with probability $q = q(w, L, \mu)$

So each read is an **informative pathogen read** with probability $p \times q$, independently of all other reads.

### The probability of detecting the pathogen

Detection happens if **at least one** read is an informative pathogen read.

The complement (missing the pathogen) means **all $T$ reads** are non-informative pathogen reads. The probability of any single read NOT being an informative pathogen read is $1 - pq$. Since reads are independent:

$$\Pr(\text{miss}) = (1 - pq)^T$$

### Solving for the minimum T

We want $\Pr(\text{miss}) \leq \delta$:

$$(1 - pq)^T \leq \delta$$

Taking the natural logarithm ($\ln$) of both sides. Remember: $\ln$ is the inverse of $e^x$. Also, $\ln$ of a number less than 1 is negative, so when we take $\ln$ of both sides, the inequality stays the same direction (both sides are negative):

$$T \cdot \ln(1 - pq) \leq \ln(\delta)$$

Now, $\ln(1-pq)$ is a **negative** number (since $pq > 0$ and we're taking $\ln$ of something less than 1). Dividing both sides by a negative number **flips the inequality**:

$$T \geq \frac{\ln(\delta)}{\ln(1-pq)} = \frac{\ln(1/\delta)}{-\ln(1-pq)}$$

Since $T$ must be a whole number (you can't do 0.7 reads), we take the **ceiling** — round up to the nearest integer:

$$T^{\ast} = \left\lceil \frac{\ln(1/\delta)}{-\ln(1-pq)} \right\rceil$$

The symbol $\lceil x \rceil$ means "ceiling of $x$" — the smallest integer greater than or equal to $x$. For example, $\lceil 2.3 \rceil = 3$ and $\lceil 4.0 \rceil = 4$.

### The approximation

When $pq$ is small (which it usually is, since $p$ is often tiny like $10^{-4}$), we can use the approximation:

$$-\ln(1-x) \approx x \quad \text{when } x \text{ is small}$$

This gives the simpler approximation:

$$T^{\ast} \approx \left\lceil \frac{\ln(1/\delta)}{pq} \right\rceil$$

This is easier to interpret: more abundance ($p$ larger) → less depth needed. Better informative probability ($q$ larger) → less depth needed. Stricter error requirement ($\delta$ smaller) → $\ln(1/\delta)$ larger → more depth needed.

### What "tight within a factor of 2" means

The exact formula and the approximation differ by at most a factor of 2 whenever $pq \leq 0.5$. This means the approximation never over- or underestimates by more than 2×, which is acceptable for a depth guideline.

### The monotonicity corollary

**Corollary:** $T^{\ast}$ increases when $\mu$ increases (more divergence = more reads needed), and increases when $w$ increases (longer window = harder to find a match = more reads needed).

**Proof:** We need to show $q(w, L, \mu)$ decreases when either $\mu$ or $w$ increases.

- More divergence ($\mu$ larger): more letters are mutated, so a run of $w$ consecutive matches is less likely. $q$ decreases. Since $T^{\ast}$ goes down as $q$ goes up, $T^{\ast}$ increases as $\mu$ increases.

- Longer window ($w$ larger): you need more consecutive letters to ALL match. This is harder for any $\mu > 0$. $q$ decreases, so $T^{\ast}$ increases.

This proves the corollary. It tells us something intuitive: harder detection situations (more diverged pathogen, or stricter matching requirement) always require more sequencing depth.

---

## 9. Theorem 2: Critical Divergence

*This is **Theorem 3.3** in the paper, Eq. 7.*

### What is $\mu_c$?

The critical divergence $\mu_c$ is the **maximum amount of divergence** a pathogen can have before it becomes impossible to guarantee detection with any practical sequencing budget $T_{\max}$.

Think of it as the boundary between Tier 2 (detectable) and Tier 3 (effectively novel).

### How is it defined?

$$\mu_c = \sup\bigl\lbrace \mu \in [0,1] : T^{\ast}(p,\mu,w,L,\delta) \leq T_{\max} \bigr\rbrace$$

Breaking this down:
- $\sup$ means "supremum" — the least upper bound, essentially the "maximum" of the set
- The set is all values of $\mu$ for which the required depth doesn't exceed our budget $T_{\max}$
- So $\mu_c$ is the largest divergence for which detection is still feasible

### How to compute it

From Theorem 1, $T^{\ast} \leq T_{\max}$ iff:

$$\frac{\ln(1/\delta)}{-\ln(1-pq)} \leq T_{\max}$$

Rearranging (the $\lceil \cdot \rceil$ only ever rounds up, so dropping it is the conservative direction):

$$-\ln(1-pq) \geq \frac{\ln(1/\delta)}{T_{\max}} \iff 1 - pq \leq \delta^{1/T_{\max}} \iff p \cdot q \geq 1 - \delta^{1/T_{\max}}$$

The critical divergence $\mu_c$ is the value of $\mu$ where that inequality becomes an equality:

$$p \cdot q(w, L, \mu_c) = 1 - \delta^{1/T_{\max}}$$

Using the approximation $1 - \delta^{1/T_{\max}} \approx \ln(1/\delta)/T_{\max}$, valid when $T_{\max}$ is large (here $1 - 0.01^{1/100000} = 4.60506 \times 10^{-5}$ against $\ln(100)/100{,}000 = 4.60517 \times 10^{-5}$, a difference of 0.002%), this is the more memorable form:

$$p \cdot q(w, L, \mu_c) \approx \frac{\ln(1/\delta)}{T_{\max}}$$

This equation says: $\mu_c$ is the divergence at which the "effective detection rate per read" $p \cdot q$ exactly equals the critical threshold $\ln(1/\delta)/T_{\max}$.

### Why does $\mu_c$ exist and is unique?

Three facts about $q(w, L, \mu)$:
1. At $\mu = 0$: $q = 1$ (no divergence, all reads are informative). So $p \cdot q = p > 0$.
2. At $\mu = 1$: $q = 0$ (all letters mutated, no reads are informative). So $p \cdot q = 0$.
3. $q$ decreases **continuously** as $\mu$ increases (no sudden jumps).

Since $p \cdot q$ starts above the threshold $\frac{\ln(1/\delta)}{T_{\max}}$ and ends below it, and moves continuously — by the **Intermediate Value Theorem** (a fundamental result in calculus: if a continuous function goes from value $A$ to value $B$, it must pass through every value in between) — there is exactly one $\mu_c$ where $p \cdot q$ equals the threshold.

### The practical formula

Using the Poisson approximation for $q$:

$$\mu_c \approx 1 - \left(\frac{-\ln\bigl(1 - \frac{\ln(1/\delta)}{p \, T_{\max}}\bigr)}{L - w + 1}\right)^{1/w}$$

This looks intimidating but is just algebra: start from $p \cdot (1 - e^{-\lambda}) = \text{threshold}$ where $\lambda = (L-w+1)(1-\mu_c)^w$, then solve for $\mu_c$.

### What $\mu_c$ tells you in practice

For our default parameters ($p = 10^{-4}$, $w = 24$, $L = 150$, $\delta = 0.01$, $T_{\max} = 100{,}000$):

$$\mu_c = 0.1349$$

This means: if a pathogen differs from the nearest reference by more than 13.49% of its genome, a 100k-read sequencing run with $w=24$ cannot guarantee detection.

**Biological interpretation:** SARS-CoV-2 vs Omicron ($\mu \approx 0.0017$) is comfortably in Tier 2. SARS-CoV-2 vs SARS-CoV-1 ($\mu \approx 0.20$) is in Tier 3. This matches what happened clinically — Omicron was rapidly detected via mNGS; a truly novel coronavirus would likely be missed.

---

## 10. Theorem 3: Optimal Window Size

*This is **Theorem 3.4** in the paper, Eq. 9.*

### The tradeoff

Choosing $w$ (the k-mer length) involves a fundamental tradeoff:

- **Smaller $w$:** Each window is shorter, so even with mutations, a read is more likely to contain an exact match. $q$ is higher → $T^{\ast}$ is lower → less sequencing needed. But shorter k-mers are also less specific — a random read from a completely unrelated organism might accidentally contain a k-mer that matches the pathogen database. More false positives.

- **Larger $w$:** Each window is longer, requiring more consecutive exact matches. Harder to satisfy with a diverged pathogen, so $q$ is lower → more sequencing needed. But longer k-mers are much more specific — they're unlikely to match by accident.

### The false-positive constraint

A **false positive** occurs when the classifier mistakenly identifies an innocent (non-pathogen) read as coming from the pathogen. This happens when a k-mer from the background DNA accidentally matches the pathogen reference.

How often does this happen? If we model the database as containing $Q \times n_{\mathrm{avg}}$ total k-mers distributed among all possible $4^w$ strings (there are $4^w$ possible strings of length $w$ over the 4-letter alphabet A, T, G, C), then the probability a random k-mer matches:

$$P_{\mathrm{FP}}(w) \approx \frac{Q \cdot n_{\mathrm{avg}}}{4^w}$$

**Why $4^w$?** Because there are 4 choices (A, T, G, C) for each of the $w$ positions, giving $4 \times 4 \times \cdots \times 4 = 4^w$ possible k-mers total.

We want this to be at most $\alpha$ (our false-positive target):

$$\frac{Q \cdot n_{\mathrm{avg}}}{4^w} \leq \alpha$$

Solving for $w$:

$$4^w \geq \frac{Q \cdot n_{\mathrm{avg}}}{\alpha}$$

Taking $\log_4$ of both sides ($\log_4(x)$ is the power you'd raise 4 to in order to get $x$):

$$w \geq \log_4\left(\frac{Q \cdot n_{\mathrm{avg}}}{\alpha}\right)$$

The smallest integer $w$ satisfying this is:

$$w^{\ast} = \left\lceil \log_4\left(\frac{Q \cdot n_{\mathrm{avg}}}{\alpha}\right) \right\rceil$$

### Why $w^{\ast}$ is optimal

We want to minimize $T^{\ast}$ (use as little sequencing as possible) while satisfying the false-positive constraint ($w \geq w^{\ast}$).

From the Corollary of Theorem 1: $T^{\ast}$ is **strictly increasing in $w$** (larger window = more sequencing needed). So among all valid $w \geq w^{\ast}$, the minimum $T^{\ast}$ is achieved at the smallest valid $w$, which is exactly $w^{\ast}$.

Any $w > w^{\ast}$ gives *even stronger* false-positive protection — $P_{\mathrm{FP}}$ keeps falling by a factor of 4 per extra base — but you have already met the target $\alpha$ at $w^{\ast}$, so that extra protection buys you nothing you asked for while costing real sequencing depth. Any $w < w^{\ast}$ would cause too many false positives.

### The sharpening tradeoff

The cost of over-specifying $w$ grows with $\mu$. For $w^{\ast} = 24$ and $\Delta = 10$, the ratio $T^{\ast}(w^{\ast}+\Delta) / T^{\ast}(w^{\ast})$ is about:

| $\mu$ | Ratio |
|---|---|
| $0$ | 1.00 (no penalty at all — $q = 1$ for every $w \leq L$) |
| $0.10$ | 2.33 |
| $\mu_c = 0.1349$ | 3.93 |
| $0.20$ | 9.71 |

So at $\mu = 0$ the choice of $w$ is free, and by $\mu = 0.20$ a window 10 bases too long costs you nearly ten times the sequencing. Note that the ratio is **finite** at $\mu_c$ — it is 3.93, not unbounded. $\mu_c$ is a statement about your depth *budget* $T_{\max}$, not a point where $q$ collapses to zero; $q(24, 150, \mu_c) = 0.4605$, which is nowhere near zero. The ratio does grow without bound eventually, but only as $\mu \to 1$, far outside the regime anyone plans an assay in.

The practical upshot: the closer you are to the detection limit, the more a too-long window costs you — but it costs you a knowable, finite multiple, which is exactly what makes it a design tradeoff rather than a cliff.

---

## 11. Real Number Walk-Through

Let's compute everything numerically for a concrete example.

### Setup

- Database: $Q = 5{,}000$ sequences, average length $n_{\mathrm{avg}} = 15{,}000$ bp
- False-positive target: $\alpha = 10^{-6}$
- Pathogen: SARS-CoV-2 Wuhan-like, compared against Omicron ($\mu = 0.0017$)
- Read length: $L = 150$ bp
- Abundance: $p = 10^{-4}$ (1 in 10,000 DNA molecules from pathogen)
- Target false-negative rate: $\delta = 0.01$ (1% chance of missing)

### Step 1: Compute $w^{\ast}$

$$w^{\ast} = \left\lceil \log_4\left(\frac{5000 \times 15000}{10^{-6}}\right) \right\rceil = \left\lceil \log_4(7.5 \times 10^{13}) \right\rceil$$

$\log_4(7.5 \times 10^{13}) = \frac{\ln(7.5 \times 10^{13})}{\ln(4)} = \frac{31.95}{1.386} = 23.05$

$$w^{\ast} = \lceil 23.05 \rceil = 24$$

### Step 2: Compute q(24, 150, 0.0017)

Using the exact recurrence (as computed by the simulator):

$$q(24, 150, 0.0017) \approx 1.00000$$

This makes sense: with $\mu = 0.0017$ (only 0.17% of letters mutated), almost every 24-letter window will be an exact match. The probability of finding no match in 127 windows is essentially zero.

### Step 3: Compute $T^{\ast}$

The exact formula uses $-\ln(1-pq)$ in the denominator:

$$T^{\ast} = \left\lceil \frac{\ln(1/0.01)}{-\ln(1 - 0.0001 \times 1.0)} \right\rceil = \left\lceil \frac{4.60517}{1.00005 \times 10^{-4}} \right\rceil = \lceil 46{,}049.7 \rceil = 46{,}050$$

It is tempting to simplify $-\ln(1-pq) \approx pq$, which gives

$$\left\lceil \frac{\ln(100)}{0.0001} \right\rceil = \lceil 46{,}051.7 \rceil = 46{,}052$$

— two reads higher. That approximation is always an *over*-estimate (since $-\ln(1-x) > x$), so it is safe to use, but the exact value is what the paper and the simulator report. The two differ by 0.004% here; the gap widens as $pq$ grows.

**Conclusion:** To guarantee detection of Omicron at 1% false-negative rate, with the given abundance and parameters, you need **46,050 reads** minimum — matching the $\mu = 0$ row of paper Table 1, since $q \approx 1$ at this divergence.

### Step 4: Compute $\mu_c$

At $T_{\max} = 100{,}000$:

$$p \cdot q(24, 150, \mu_c) = \frac{\ln(100)}{100{,}000} = \frac{4.605}{100{,}000} = 0.00004605$$

$$q(24, 150, \mu_c) = \frac{0.00004605}{0.0001} = 0.4605$$

We need to find $\mu_c$ such that the recurrence gives $P_{150}(24, 1-\mu_c) = 1 - 0.4605 = 0.5395$.

Using numerical bisection (as implemented in the simulator): $\mu_c \approx 0.1349$.

---

## 12. Why the Poisson Approximation Fails

### What is the Poisson approximation?

The **Poisson approximation** is a shortcut formula for $q$:

$$q(w, L, \mu) \approx 1 - e^{-(L-w+1)(1-\mu)^w}$$

Here $e \approx 2.718$ is Euler's number, and $(L-w+1)(1-\mu)^w$ is the expected number of matching windows (since there are $L-w+1$ windows, each matching with probability $(1-\mu)^w$).

This approximation comes from treating the $L-w+1$ windows as **independent**. If you have many independent events each with a small probability, the total number follows a Poisson distribution.

### When does it fail?

The problem is that windows are **not independent**. Window 1 (positions 1–24) and Window 2 (positions 2–25) share 23 out of 24 letters. If window 1 is a match, window 2 is almost certainly also a match. They're highly correlated.

It is tempting to think the approximation is safe whenever $(1-\mu)^w$ is small, so that matching windows are rare. That is not what rescues it. The matching windows **clump**: a run of $m \ge w$ consecutive matching bases produces $m-w+1$ matching windows all at once, so the expected count $\lambda = (L-w+1)(1-\mu)^w$ counts one lucky run many times over — by a factor of roughly $1/\mu$. Shrinking $\lambda$ does not fix that. At $\mu = 0.30$, $\lambda \approx 0.024$ is tiny, and the approximation still overestimates $q$ by a factor of 3.2.

What actually controls the error is $\mu$ itself. The overestimate is harmless at low divergence and gets worse as divergence rises:

| $\mu$ | exact $q$ | Poisson $q$ | ratio |
|---|---|---|---|
| 0.05 | 0.986 | 1.000 | 1.01 |
| 0.10 | 0.744 | 1.000 | 1.34 |
| 0.12 | 0.580 | 0.997 | 1.72 |
| 0.20 | 0.119 | 0.451 | 3.80 |

So the danger zone is not the low-$\mu$ regime — there the approximation is fine — but the moderate-to-high divergence regime at and above $\mu \approx 0.10$, which is exactly where the detection guarantee is marginal. A lab using the Poisson formula at $\mu = 0.12$ would plan for 46,175 reads when 79,409 are needed, and would compute a critical divergence of 0.199 instead of the true 0.1349, treating a $\mu = 0.18$ pathogen as comfortably detectable when it actually needs 243,553 reads. That is a false sense of security precisely where it matters most.

**This is why the exact linear recurrence is essential** — it correctly handles the correlation between overlapping windows.

---

## 13. Glossary

Every technical term used in this project, defined simply.

**Abundance ($p$):** The fraction of all DNA molecules in a sample that come from a specific pathogen. If 1 in every 10,000 molecules is from the pathogen, $p = 0.0001$.

**Base pair (bp):** One letter in the DNA alphabet (A, T, G, or C). Genome lengths are measured in base pairs.

**Bernoulli trial:** A single random experiment with exactly two outcomes: success (with probability $p$) and failure (with probability $1-p$). Named after mathematician Jacob Bernoulli.

**Ceiling function ($\lceil x \rceil$):** Rounds a number up to the nearest integer. $\lceil 2.1 \rceil = 3$, $\lceil 5.0 \rceil = 5$.

**Classifier:** An algorithm that takes a DNA read and assigns it to a species or pathogen. K-mer classifiers work by looking for exact substring matches in a reference database.

**Corollary:** A result that follows easily from a theorem. If the theorem is the main discovery, the corollary is a direct consequence that requires little additional work to prove.

**Critical divergence ($\mu_c$):** The maximum divergence for which detection can be guaranteed with a given sequencing budget. Pathogens with $\mu > \mu_c$ cannot be reliably caught using reference-based methods.

**Divergence ($\mu$):** The per-base rate of mutation between a pathogen genome and the nearest reference. $\mu = 0.05$ means 5% of letters have changed.

**False negative:** Concluding that a pathogen is absent when it is actually present. The false-negative rate $\delta$ is the acceptable probability of this error.

**False positive:** Concluding that a pathogen is present when it is actually absent. Controlled by the parameter $\alpha$.

**Floor function ($\lfloor x \rfloor$):** Rounds down to nearest integer. The opposite of ceiling.

**Genome:** The complete set of DNA instructions for an organism, from viruses (a few thousand letters) to humans (3 billion letters).

**Informative read:** A sequencing read that contains at least one k-mer exactly matching the reference database, allowing a classifier to recognize it.

**Intermediate Value Theorem (IVT):** If a continuous function takes value $A$ at one point and value $B$ at another, it must take every value between $A$ and $B$ somewhere in between. Fundamental result from calculus.

**k-mer (or w-mer):** A substring of length $k$ (or $w$). "ATGCGA" contains the 3-mers ATG, TGC, GCG, CGA.

**$\ln$ (natural logarithm):** The inverse of $e^x$. $\ln(e) = 1$, $\ln(1) = 0$, $\ln(x) < 0$ for $x < 1$. Used extensively because exponential decay ($e^{-\lambda}$) appears naturally in probability.

**$\log_4$:** Logarithm base 4. $\log_4(16) = 2$ because $4^2 = 16$. Appears because the DNA alphabet has 4 letters.

**mNGS (metagenomic next-generation sequencing):** Sequencing all DNA in a sample to identify organisms without prior knowledge of what's there.

**Monotone / monotonically increasing:** A function is monotonically increasing if it always goes up as the input increases. $T^{\ast}$ is monotonically increasing in $\mu$ (more divergence always requires more depth).

**Mutation:** A change in a DNA letter. Pathogens mutate as they replicate, causing divergence from the reference over time.

**Pathogen:** A disease-causing organism — bacteria, viruses, fungi, parasites.

**Poisson approximation:** A shortcut that approximates the probability of at least one event occurring from many nearly-independent, low-probability trials. Valid when the events are roughly independent; breaks down when windows overlap heavily.

**Read:** A short DNA sequence produced by a sequencing machine, typically 75–300 base pairs long.

**Recurrence relation:** A formula that defines each term of a sequence using previous terms. $P_n = f(P_{n-1}, P_{n-2}, \ldots)$. Efficient to compute because you build up the answer step by step.

**Reference database:** A collection of known genome sequences used for comparison. If a read matches something in the database, it's classified as that organism.

**Run:** A sequence of consecutive identical outcomes. A "run of $w$ successes" means $w$ consecutive Bernoulli trials all resulted in success.

**Sequencing depth ($T$):** The total number of reads produced by a sequencing run. Higher depth = more chances to observe a rare pathogen.

**$\sup$ (supremum):** The least upper bound of a set. For a set of numbers, the supremum is the smallest number that is $\geq$ every element of the set. For a finite set, it's just the maximum.

**Theorem:** A mathematical statement that has been proven to be true from first principles, using only logic and previously established results.

**Tier:** Our classification of pathogens by divergence: Tier 1 (exact match), Tier 2 (divergent but detectable), Tier 3 (effectively novel).

**$T^{\ast}$:** The minimum sequencing depth that guarantees $\delta$-detection. The central result of this paper.

**Window:** A substring of a read of length $w$, extracted by sliding along the read one position at a time. A read of length $L$ has $L - w + 1$ windows.

**$w^{\ast}$:** The optimal k-mer window length that minimizes required sequencing depth while satisfying the false-positive constraint.

---

*For the formal mathematical statements and proofs in standard notation, see [README.md](README.md). For the code that validates all results, see [code/simulator.py](code/simulator.py).*