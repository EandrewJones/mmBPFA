# Learning the number of factors: a new loadings prior for mmBPFA

Status: draft, 2026-10-19. Two choices are still waiting on running experiments: candidate C versus E3, and the spike width. See "Open questions".

## 1. Intent

mmBPFA is a Bayesian Gaussian-copula factor model. It uses Hoff's (2007) extended rank likelihood (ERL), so it can scale any mix of binary, ordinal, count and continuous items without a model for their margins. Two properties are meant to set it apart, and both are contributions the author relies on:

1. **The number of dimensions K is learned from the data.** A user who has no idea what K is should get an answer they can trust: one that comes from fit to the data, not from the starting value or from hidden prior choices.
2. **The fit says which items load on which dimension.** That item-level map is how dimensions get meaning ("these votes and these words load here"). Without it, the dimensions can't be interpreted.

The method is meant to be general: multi-modal item-response models across many applications, not only the UN General Assembly paper. A later goal is to scale to very large N×P with sparse observations, such as an online recommender over multi-modal rankings, using variational or stochastic variational inference. This design must not close that door.

### Why this change is needed

The current two-parameter IBP prior (Knowles & Ghahramani 2011) with slab λ_jk ~ N(0, 1/γ_k), γ_k ~ Gamma(1, rate d = 100) does not learn K reliably once the earlier sampler bugs are fixed:

- **The slab is badly calibrated on the identified scale.** The ERL identifies only the copula correlation C = D^{-1/2}(ΛΛ′ + I)D^{-1/2}, i.e. copula loadings r_j = λ_j / √(1 + ‖λ_j‖²). With d = 100, 86% of the prior mass for an active loading has |r| > 0.9. Binary probit data can barely tell λ = 2 from λ = 20, so spurious "doublet" factors settle at large loadings and persist.
- **The sampler cannot prune.** The posterior does not depend on K_start, yet runs from different starts disagree badly (for example K = 2 versus K = 16–32 on the same data). Columns can only empty one element at a time against rich-get-richer prior odds, and new features are born only as singletons.
- **Small factors are not identified.** One-item factors are exactly unidentified under the ERL. Two-item factors are weakly identified: they are confounded with a residual correlation.
- **Theory gives no protection.** IBP-type priors are consistent for K only when the true sparsity is known (Ohn, Lin & Kim 2023), and under the IBP the expected K grows with the number of items (Ohn 2022).

The literature review is in the session scratchpad (`priors-literature.md`); its main points are cited inline here.

## 2. Goals, non-goals, anti-goals

**Goals**
- K estimates that agree across chains started from widely different K (for example 2, 10, 30) and that recover the true K in simulations whose factors are detectable.
- An item-level inclusion map z_jk (which item loads on which dimension), recovered in simulations.
- Priors stated and calibrated on the copula-correlation scale that users interpret.
- Copula-correlation recovery at least as good as the current fixed sampler.
- A model structure that maps onto closed-form variational inference later: finite truncation, conjugate-exponential blocks, and no correctness that depends on collapsed-only moves.

**Non-goals (for this change)**
- Implementing VI or SVI. Only keep the structure compatible.
- A C++ performance rewrite of the sweep. That is a separate, planned piece of work. New code should be easy to port: plain linear algebra, no R-only tricks inside the inner loops.
- A dynamic (over time) model.
- Changing the ERL. Inference under the ERL for VI is an open problem and is not addressed here.

**Anti-goals: things that would make this a failure even if it technically runs**
- A sampler that reports a confident K which depends on where it started. This is today's failure, and it is worse than being slow.
- Losing the item-level map. For example, a column-only shrinkage prior such as plain CUSP.
- Hiding the influence of priors. Every prior that affects K must be visible, calibrated on an interpretable scale, and documented.
- Counting unidentified factors (one- or two-item) as dimensions.

## 3. Decisions already made

These were decided by the author.
- The item-level map is essential, so the prior must have element-wise inclusion indicators.
- Learning K from data is a core contribution. A reliable fit to data matters more than keeping the IBP for its own sake.
- Column-level cumulative shrinkage (CUSP, Legramanti, Durante & Dunson 2020) is the right family to build on.
- Try the Kowal & Canale (2023) parameter expansion to improve mixing over K.
- Keep the method general (any mix of margins, large and sparse data), not tuned to the UN paper.
- **The acceptance criteria in §7 are soft.** They guide judgement, and any miss is reported with its size and likely cause, not treated as a blocker.
- **The IBP is removed in this change.** It performed worse in the spike: even calibrated, it was the only candidate whose chains disagreed across starts. The package is private, so there is no deprecation period. Tagged releases preserve it (see §6.1).
- **Users must be able to retrieve every prior setting.** This includes defaults and anything derived from a calibration (see §6.1).

## 4. The model

Unchanged parts:
- the latent data x_ij = Σ_k λ_jk ω_ik − α_j + ε_ij, with ε ~ N(0, 1);
- ω_i ~ N(0, I);
- a flat prior on α_j;
- the ERL on x.

Truncation H is the maximum number of columns (see §6.4). Columns are ordered, and k indexes them.

### 4.1 Column level: cumulative shrinkage (shared by C and E3)

- Each column has a categorical indicator c_k ∈ {1, …, H}. Column k is **on** (slab) when c_k > k and **off** (spike) otherwise.
- P(c_k = l) = w_l, with stick-breaking weights w_l = v_l Π_{m<l}(1 − v_m), v_l ~ Beta(1, a_cusp) for l < H, and v_H = 1.
- This makes the prior probability of a column being off increase with its index: π_k = Σ_{l≤k} w_l. It is exactly the IBP stick-breaking weight applied per column (Frühwirth-Schnatter 2023).
- **a_cusp ~ Gamma(1, 1) is learned**, so the expected number of active columns is not fixed by the user. Its full conditional is conjugate: Gamma(1 + H − 1, 1 − Σ_{l<H} log(1 − v_l)).

### 4.2 Item level: inclusion within columns (shared by C and E3)

- z_jk ~ Bern(ρ_k), with ρ_k ~ Beta(1, 1) per column.
- Only z_jk = 1 entries carry a loading.
- This is the Schiavon, Canale & Dunson (2022) structured increasing shrinkage (SIS) construction without meta-covariates. Their default scales the inclusion probability by c_p = 2e·log(p)/p. That is left as an open question (§9).

### 4.3 Loadings: two candidates

**Candidate C (SIS-style, no expansion)**
- λ_jk | z_jk = 1 ~ N(0, θ_k).
- θ_k = θ∞ when the column is off. When it is on, θ_k ~ InvGamma(1, ½).
- **Calibration.** For a single active factor, InvGamma(1, ½) makes the copula loading r = λ/√(1 + λ²) exactly Uniform(−1, 1). This follows because λ ~ t₂(0, 1/√2) and r² ~ Beta(½, 1); it was verified by simulation and independently derived.
- Spike: θ∞ = 0.05 in the spike runs (the CUSP default, used here as a variance).
- **Column-indicator step:** conditional on the active loadings of column k. The spike likelihood is Π N(λ_jk; 0, θ∞). The slab likelihood is the multivariate t from integrating θ_k out of the InvGamma(1, ½) slab.

**Candidate E3 (Kowal–Canale expansion with an absolute spike)**
- λ_jk = z_jk · η_k · ξ_jk.
- Column scale: η_k ~ N(0, s_k²) when the column is on, and η_k ~ N(0, σ_spike²) when off. **The spike is absolute**, not relative to s_k².
- s_k⁻² ~ Gamma(1, rate 0.35). This is calibrated so that, for a single active factor, P(|r| > 0.9) ≈ 0.10, matching the uniform tail. The median of |r| is about 0.4 rather than 0.5.
- Item multipliers: ξ_jk ~ N(m_jk, 1), with m_jk ∈ {−1, +1} equally likely.
- **Column-indicator step:** depends only on the scalar η_k, comparing N(η_k; 0, σ_spike²) against N(η_k; 0, s_k²). That is why it mixes.
- Two literal variants were tried and rejected (§8):
  - Kowal–Canale's per-column relative spike: the evidence for a column being on is capped.
  - A global scale: K = 0 becomes absorbing.

### 4.4 What counts as a dimension

- K for a draw is the number of **on** columns with **at least 3 active items**.
- This is a necessary part of the identification condition, not a heuristic. Factor identification needs any q columns to have at least 2q + 1 non-zero rows, so each column needs at least 3 (the "3579 rule", Frühwirth-Schnatter, Hosszejni & Lopes). It also implies K ≤ (P − 1)/2.
- Builder's choice: also report the share of draws that pass the full 3579 rule, as a diagnostic.
- Report the posterior distribution of K, not only the median.

### 4.5 What gets reported (my call)

- **The model's correlation uses all H columns.** The copula correlation, the likelihood and GMP are computed from every column, including off columns and on columns with fewer than 3 items, because that is the model's posterior. Criterion 4 (§7) compares this full-model correlation.
- **Loadings and dimensions use only the counted columns.** Users see loadings, ideal points and item maps for the counted dimensions of §4.4 only.
- **The gap is reported.** Each summary also reports the share of the correlation that comes from uncounted columns, the mean absolute difference between the full-model C and the C rebuilt from counted columns alone. That makes it visible when uncounted columns matter.
- **Storage.** Stored draws keep all H columns (P×H per draw, plus the on/off state) so both views can be computed afterwards. Relabelling (`align_dimensions`, `find_stable_dims`) works on the counted columns.

## 5. Gibbs sweep

The order within a sweep is fixed. Everything not listed is unchanged from the current fixed sampler.

1. **x:** the ERL truncated-normal step, sequential by level (already fixed).
2. **Loadings sweep, column by column.** Use running residuals: add column k back, update it, subtract it.
   - **C:** for each item, a collapsed Bernoulli for z_jk with λ_jk integrated out under N(0, θ_k), then λ_jk | z = 1 from its Gaussian conditional.
   - **E3:**
     - z_jk with λ_jk | z = 1 ~ N(η_k m_jk, η_k²) integrated out. The log-odds are log ρ_k/(1 − ρ_k) − ½ log(η_k² p_k) + ½ p_k μ_jk² − ½ m_jk², where p_k = ω_kᵀω_k + η_k⁻² and μ_jk = (ω_kᵀe_j + m_jk/η_k)/p_k.
     - Then λ_jk, then ξ_jk = λ_jk/η_k, drawing from the prior when z = 0.
     - Then m_jk with P(m = +1) = logit⁻¹(2ξ_jk).
     - Then η_k from its Gaussian conditional given ξ and z.
3. **ρ_k:** Beta(1 + m_k, 1 + P − m_k).
4. **Column indicators c_k:** as in §4.3, then v_l, then a_cusp.
5. **Column variances:**
   - C: θ_k from its InvGamma conditional when the column is on, θ∞ when off.
   - E3: s_k² from its InvGamma conditional when on, and from its prior when off.
6. **ω** (existing step, now with no standardization), then **α** (existing, now exact).

The existing IBP pieces are removed, along with their tests:
- the IBP z step and singleton new-feature move;
- the IBP a and b updates;
- `sample_gamma_k`;
- `remove_weakly_active`, which is replaced by §4.4.

The pre-sparse phase becomes unnecessary: CUSP needs no warm start without sparsity. The builder should confirm this in the recovery tests before deleting it; otherwise keep it as an option.

## 6. Constraints and interfaces

### 6.1 Package interface
- `mmBPFA_sampler()` keeps its signature where it can.
- `K_start` now means the number of columns that start **on**.
- Add an argument for the truncation H, with the default in §6.4.
- Replace `share_noise` with prior arguments for the slab and spike. How users set priors is my call: a `priors` argument that takes either a named preset or a list of raw hyperparameters. The default preset is `"uniform"`, giving a uniform prior on single-factor copula loadings. Raw values override the preset.
- **Every prior is retrievable** (author's requirement). The sampler output carries a complete, flat record of every hyperparameter actually used: the slab and spike; the column inclusion prior; the a_cusp hyperprior; truncation H; and anything derived from a preset, with the preset name. `get_priors(results)` returns it for both raw and processed results, and print/summary methods show it.
- **Old samplers live on as tagged releases** (the tag names are my call):
  - `v1.0-paper` tags the original master commit (761a9ba) that produced the published results, bugs included;
  - `v1.1-ibp` tags the corrected IBP sampler once the `fix/sampler-correctness` branch is merged.
  The README / NEWS points to both. Creating and pushing the tags needs the author's go-ahead at the time.
- Stored draws must include, per draw: the on/off state of each column, z, λ (on the sampler scale), ω, α, and the cutpoints. Downstream code must not depend on the IBP-only fields.

### 6.2 Post-processing
- `process_results()` keeps pivot-based relabelling (`align_dimensions`), since CUSP ordering reduces label switching but does not remove it.
- "Copula Loadings" stays.
- The K summary uses §4.4.

### 6.3 Variational compatibility
- Every conditional above is conjugate given the finite truncation H. C maps onto Legramanti's (2020) CUSP variational Bayes; E3 adds Gaussian q(η) and q(ξ) factors and a categorical q(m).
- No step may rely on a collapsed move for correctness.

### 6.4 Truncation
- Default H = min(⌊(P − 1)/2⌋, 50).
- If `K_start` exceeds H, cap it at H with a message (my call).
- If a chain spends more than a small share of draws with every column on, warn the user to raise H.
- CUSP's adaptive truncation (exp(−1 − 5e-4·t)) is an allowed alternative. It is the builder's call if fixed H proves too slow at large P.

### 6.5 Engineering
- R and RcppArmadillo, as now.
- Every sampler step is developed test-first with statistical tests, following the repository's conventions in `CLAUDE.md`:
  - conditional samplers checked against analytic posteriors using `expect_z` and `expect_rel`;
  - prior-invariance runs with no data;
  - recovery runs on data simulated from the model.

## 7. Validation: acceptance criteria

These criteria are soft (author's decision). They guide the choice, and any miss is reported with its size and likely cause. The targets, on the simulation grid below (two data seeds per setting, three starts each, 2,000 iterations with 1,000 burn-in):

1. **Agreement:** the median K from starts 2, 10 and 30 differ by at most 1 in every setting.
2. **Accuracy:** the median K is within 1 of the true K in every setting, and equal to it in at least 3/4 of them.
3. **Item map:** F1 for item inclusion is at least that of the current fixed sampler in each setting.
4. **Correlation:** copula-correlation error is no worse than the current fixed sampler's, within 0.005.
5. **Prior invariance:** with no data, the sampler reproduces the prior distribution of the number of on columns.

Grid: N = 200; P ∈ {60, 200}; K ∈ {2, 5}; binary-only and binary-plus-count margins; 5% missing at random.

**Truncation in the grid (my call).** The acceptance runs set H = 30 explicitly, overriding the default, which is 29 at P = 60. "Start at 30" then means "every column on", the most adverse start available.

**Baseline for criteria 3 and 4 (my call).** The baseline is spike candidate A: the collapsed IBP with the calibrated slab, run through the same harness, settings and data. It is the strongest IBP comparator, and the package defaults (d = 100) are already known to fail. Because the IBP is removed in this change, its grid results are saved as a fixed reference file next to the ported recovery test, together with the `v1.1-ibp` commit that produced them.

**Mixing stage (my call).** Medians alone can't tell real mixing from churn among weak columns, so the C-versus-E3 decision also uses:
- the effective sample size of K and of the leading copula correlations per second of runtime, on the grid;
- one larger run at N = 200, P = 1000 (mixed margins, K = 5, three starts, 5,000 iterations), to check that the mixing advantage holds as P grows. The harness lives in the scratchpad spike (`spike/harness.R`) and should be ported into `tests/` as a slow, `skip_on_cran` recovery test.

Unit-level tests (each written first and watched failing):
- each conditional against its analytic posterior;
- the column-indicator step against a hand-computed categorical;
- the calibration claim, P(|r| > 0.9) under the slab prior;
- the K definition (on columns with at least 3 items).

## 7b. Stage 2: stress tests on diverse data-generating processes

Run once the field is down to one or two candidates (requested by the author). The point is to check the "general method" claim against the kinds of data that arrive in ranking and recommender systems, not only in roll-call votes.

### Data-generating processes

Each process is simulated with known latent structure where that structure is defined.

| Family | Process | What it stresses |
|---|---|---|
| Binary | Balanced votes; **extreme imbalance** (click-through rates of 1–5%) | Thresholds deep in the tails; weak information per cell |
| Ordinal | Likert 1–5 with item-specific thresholds; **J-shaped** ratings (mostly 1s and 5s); response-style shifts per user | Heavy ties; a user effect that looks like a factor |
| Nominal | **Unordered multinomial choices** (pick 1 of C), one-hot encoded | One-hot columns are negatively dependent by construction and can create spurious factors |
| Counts | Poisson; negative binomial (overdispersed); **zero-inflated** (most cells 0); power-law popularity | Massive ties at 0; heavy tails |
| Continuous | Normal; lognormal (dwell time); Student-t (outliers); Beta (completion rates); **censored** (watch time capped) | Margin-free handling of skew, bounds and censoring |
| Rankings | **Top-k lists and pairwise comparisons** from Plackett–Luce / Thurstone models | Orderings *within a user (row)*, unlike the ERL's within-item (column) ordering |
| Multi-view | The same items observed in several modalities (rating + click + text count), with shared factors and **modality-specific "method" factors** | Whether modality-only dimensions are found and kept apart |
| Missingness | MCAR at 50–99.9% sparsity; **missing not at random** (exposure bias: users see popular items, and high ratings are more likely observed); power-law activity per user and per item | The ERL treats missing cells as uninformative; real recommender data violate that |
| Misspecification | Clustered or skewed users (a mixture of Gaussians for ω); t-copula tail dependence; near-duplicate items (a genuine 2-item block); random-responder rows (bots) | Whether non-Gaussian dependence inflates K; how doublets and outliers are treated |

### Two structural limits these tests will probe

1. **The ERL orders values within each item (column).** Unordered categories have no order, so they need an encoding (one-hot) or a multinomial-probit extension. Within-user rankings constrain rows, not columns, so they need either an encoding (rank position per item, NA when unranked) or a row-wise rank likelihood. The tests will show how far encodings get before an extension is needed.
2. **Missing-not-at-random data.** The model assumes missing cells carry no information. Exposure-biased data will show how badly K and the loadings drift when that fails.

Both are findings to report, not things this change has to fix.

### Metrics

- **When the latent structure is defined:** K recovery, item-map F1, and copula-correlation error, as in §7.
- **Always:**
  - agreement across starts;
  - held-out predictive performance on masked cells (log predictive probability; AUC for binary cells; NDCG@k per user for ranking-type data);
  - runtime, and effective sample size per second.
- **When the true copula is not Gaussian** (t-copula, Plackett–Luce): report rank-correlation recovery (Kendall's τ) and predictive metrics in place of C error. Report whether K inflates relative to the generating structure.
- **Reference points:** a simple baseline is reported next to each process, so improvements are interpretable. Candidates are normal-score PCA with parallel analysis for K, polychoric exploratory factor analysis, and probabilistic matrix factorization for predictive metrics. Which baselines exactly is the builder's call.

### Scale

N and P up to the low thousands, at the sparsity levels above. Larger runs wait for the C++ sweep.

## 8. Evidence so far

All runs below are throwaway spike runs on simulated data, every candidate with the calibrated slab. Each cell gives the median K from starts at 2, 10 and 30.

| Setting | A: IBP | B: IBP + feature moves | C: CUSP + item inclusion | D: finite beta–Bernoulli |
|---|---|---|---|---|
| K=2, P=60 binary | 2/3/19 | 3/3/3 | 2/2/2 | 3/3/2 |
| K=2, P=200 mixed | 2/3/2 | 3/2/3 | 2/2/3 | 3/3/2 |
| K=5, P=60 binary | 6/7/17 | 7/6/7 | 6/5/6 | 8/5/7 |
| K=5, P=200 mixed | 6/5/6 | 5.5/6/7 | 5/5/6 | 6/5/8 |

- C was the most consistent candidate.
- Correlation error for C was slightly worse (for example 0.028 against 0.020). Its relative spike leaves the off columns non-negligible.
- B cost 4–6× the runtime with no gain.

Kowal–Canale expansion variants (K = 2, P = 60, binary, 1,000 iterations):

| Variant | Median K from 2 / 10 / 30 | Share of iterations K changes |
|---|---|---|
| C | 2–3 / – / 3 | 20–27% |
| E, as published (relative spike, per-column scale) | 0 / – / 0 | |
| E2, global scale | 3 / – / 0 (absorbing) | |
| E3, absolute spike sd 0.05 | 3 / 4 / 4 | 53–60% |

### C versus E3 on the full grid

Two data seeds per setting, starts at 2, 10 and 30, 2,000 iterations.

| Metric | C (spike variance 0.05) | E3 (spike sd 0.05) |
|---|---|---|
| Starts agree within 1 (criterion 1) | 100% | 75% |
| Median K exactly right (criterion 2: ≥ 75%) | 81% | 21% |
| Median K within 1 of truth | 100% | 79% |
| Mean over-count | +0.19 | +0.94 |
| Item-map F1 | 0.80 | 0.70 |
| Copula-correlation error | 0.031 | 0.029 |
| Share of iterations where K changes | 0.25 | 0.56 |

- E3's higher rate of K changes is churn among weak columns, not useful mixing: it comes with over-counting and a worse item map.
- **Against baseline A on identical data:** C's item-map F1 is equal or better in all 8 settings (criterion 3 passes). Its correlation error is within 0.005 in 6 of 8 settings. It misses criterion 4 only at K = 2, P = 60, by 0.006–0.007, which is attributed to its non-negligible spike (sd ≈ 0.22).
**Spike-width ablations.**

| Candidate | Starts within 1 | Median K exact | F1 | Correlation error within 0.005 of baseline | K-change rate |
|---|---|---|---|---|---|
| C, spike variance 0.05 | 100% | 81% | 0.80 | 6/8 | 0.25 |
| C, spike sd 0.05 | 44% | 65% | 0.78 | 7/8 | 0.08 |
| C, spike sd 0.02 | 31% | 54% | 0.74 | 4/8 | 0.07 |
| E3, spike sd 0.05 / 0.02 / 0.01 | 56–75% | 21–29% | 0.70–0.72 | 8/8 | 0.29–0.56 |

- A narrow spike freezes K, the Kowal–Canale "stuck K" failure: columns with real loadings can't move into it.
- **Decision (Task 1):** C with spike variance 0.05.

**Caveat found during implementation.** The spike harness drew the CUSP concentration given the sticks. Once every column is off, the sticks round to 1 and the concentration to exactly 0, which makes K = 0 absorbing in floating point. The package draws the concentration with the sticks integrated out instead. So the K = 0 collapses of E and E2 may come partly from this trap, not only from their design flaws as diagnosed in the table above. Candidate C never reached K = 0, so the decision stands.

## 9. Open questions

**Waiting on running experiments** (two seeds, the full grid):
1. **C or E3.** On the full grid, C passes criteria 1–3 and E3 does not (see §8). C is the leading candidate, pending the smaller-spike runs.
2. **E3 spike width.** sd 0.05, 0.02 or 0.01. A smaller spike should cut over-counting and correlation noise, but may slow mixing.
3. **Why there's an over-count of about 1 across candidates.** Likely weak spurious factors that just meet the 3-item bar. Check whether a sharper spike or the SIS c_p scaling of inclusion probabilities removes it.

**For the author:**
4. **Naming.** C and E3 keep only the IBP's stick-breaking weights, and they apply them per column rather than per item–dimension cell. Item inclusion is a separate beta–Bernoulli layer. The model is no longer a beta process, so "beta-process factor analysis" needs rewording in the package description and the paper. This doesn't block implementation.

**Left to the builder:**
- fixed versus adaptive truncation (§6.4);
- whether to keep the pre-sparse phase;
- the exact form of diagnostics and warnings;
- how to structure the C++ for the new loadings sweep.

## 10. References

- Frühwirth-Schnatter, Hosszejni & Lopes (2024), Sparse Bayesian factor analysis when the number of factors is unknown, Bayesian Analysis; arXiv:2301.06459.
- Frühwirth-Schnatter (2023), Generalized cumulative shrinkage process priors, Phil. Trans. R. Soc. A; arXiv:2303.00473.
- Hoff (2007), Extending the rank likelihood for semiparametric copula estimation, Annals of Applied Statistics.
- Knowles & Ghahramani (2011), Nonparametric Bayesian sparse factor models, Annals of Applied Statistics.
- Kowal & Canale (2023), Semiparametric functional factor models with Bayesian rank selection, Bayesian Analysis; arXiv:2108.02151.
- Legramanti, Durante & Dunson (2020), Bayesian cumulative shrinkage for infinite factorizations, Biometrika; arXiv:1902.04349.
- Legramanti (2020), Variational Bayes for Gaussian factor models under the cumulative shrinkage process; arXiv:2008.05310.
- Murray, Dunson, Carin & Lucas (2013), Bayesian Gaussian copula factor models for mixed data, JASA.
- Ohn (2022), arXiv:2206.08002. Ohn, Lin & Kim (2023), arXiv:2305.18488.
- Schiavon, Canale & Dunson (2022), Generalized infinite factorization models, Biometrika; arXiv:2103.10333.
