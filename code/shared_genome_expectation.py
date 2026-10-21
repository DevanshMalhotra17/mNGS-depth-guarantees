import numpy as np, math
rng = np.random.default_rng(7)
w, L, p, T = 24, 150, 1e-4, 100_000

def q_exact(mu):
    s = 1 - mu
    P = [1.0] * (L + 1)
    for n in range(w, L + 1):
        P[n] = sum((s ** j) * mu * P[n - j - 1] for j in range(w))
    return 1 - P[L]

def f_of_genome(N, mu):
    # fraction of read start positions (0..N-L) whose read has >=1 exactly matching w-window
    mism = rng.random(N) < mu
    c = np.concatenate(([0], np.cumsum(mism)))
    win_ok = (c[w:] - c[:-w]) == 0                 # window j matches: j=0..N-w
    cw = np.concatenate(([0], np.cumsum(win_ok)))
    pos = np.arange(N - L + 1)
    inf = (cw[pos + (L - w) + 1] - cw[pos]) > 0
    return inf.mean()

def run(N, mu, G):
    f = np.array([f_of_genome(N, mu) for _ in range(G)])
    D_shared = np.mean(1 - (1 - p * f) ** T)       # exact E over genomes; reads sampled analytically
    se = np.std(1 - (1 - p * f) ** T) / math.sqrt(G)
    q = q_exact(mu)
    D_theory = 1 - (1 - p * q) ** T
    return D_theory, D_shared, se, f.mean(), f.std()

print("N      mu    D_theory  E[D | shared genome]  (MC se)   mean f   sd f   gap")
for N, G in ((5000, 20000), (29903, 5000)):
    for mu in (0.15, 0.18, 0.20):
        Dt, Ds, se, fm, fs = run(N, mu, G)
        print(f"{N:<6} {mu:.2f}  {Dt:.4f}    {Ds:.4f}            ({se:.4f})  {fm:.4f}  {fs:.4f}  {Dt-Ds:.4f}")
