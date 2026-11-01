# CAJD — Composite Approximate Joint Diagonalization

MATLAB code for **Composite Approximate Joint Diagonalization (CAJD)**, a blind source separation method for data in which some sources follow a **linear** model and others a **bilinear** (spatio-temporal) model. It also includes the two baselines CAJD is compared against, plain **AJD** and **bilinear AJD (BAJD)**. All three use the same Jacobi-like solver, **Gauss Planar Transformations (GPT)**.

The main application is extracting **Event-Related Potentials (ERPs)** from EEG. An ERP has a fixed spatial *and* temporal pattern (bilinear), and it sits on top of ongoing brain activity that only has a fixed spatial pattern (linear). CAJD models both together.

This is the code behind:

> L. Korczowski, F. Bouchard, C. Jutten, M. Congedo, **"Mining the Bilinear Structure of Data with Approximate Joint Diagonalization"**, *EUSIPCO 2016 – 24th European Signal Processing Conference*, Budapest, Aug. 2016, pp. 667–671. [hal-01357245](https://hal.science/hal-01357245v1) · [PDF](https://hal.science/hal-01357245/document)

![Sources estimated by CAJD on EEG: temporal patterns (left) and scalp maps (right). c06 is the P3b, c11 the early visual complex.](pic/EEG_independent_sources.png)

---

## The three models

All matrices are real. `N` is the number of sensors (and of sources), `T` the number of time samples.

| Method | Data model | Input to diagonalize | Unmixing found | Cost minimized |
|---|---|---|---|---|
| **AJD** | linear: `x(t) = A s(t) + n(t)` | symmetric target matrices `R_l = A Λ_l Aᵀ` (covariances, co-spectra…) | spatial `B` | `f(B) = Σ_l ‖off(Bᵀ R_l B)‖²_F` |
| **BAJD** | bilinear: `X_k = A S_k Eᵀ + N_k` | raw trials `X_k` (N×T), not their covariances | spatial `B` and temporal `D` | `f_b(B,D) = Σ_k ‖off(Bᵀ X_k D)‖²_F` |
| **CAJD** | both at once, sharing the same `A` | trials `X_k` **and** target matrices `R_l` | `B` and `D` | `f_c = α f(B) + (1−α) f_b(B,D)` |

- `A` (N×N) is the spatial mixing matrix and `E` (T×N) the temporal mixing matrix. Neither is assumed orthogonal.
- `S_k` and `Λ_l` are diagonal and hold the source amplitudes and powers.
- `off(·)` sets the diagonal to zero.
- The estimated sources are `Ŝ_k = Bᵀ X_k D`. Spatial patterns are `Â = B⁻¹` and temporal patterns `Ê = D⁺` (pseudo-inverse).

**Why combine them?** In real ERP data the evoked part is bilinear and the background activity is linear. Fitting only one of the two models misses the other part. The linear targets add *diversity* that helps pin down `B`, so CAJD estimates the spatial unmixing more accurately than BAJD on noisy data. For the temporal unmixing `D` the two perform about the same.

### Solver: Gauss Planar Transformations (GPT)

Starting from `B0` and `D0`, one **sweep** goes over every ordered pair `(i, j)` with `i ≠ j` and applies a non-orthogonal planar update:

```
b_i ← b_i + β b_j        d_i ← d_i + γ d_j
```

`β` and `γ` are the closed-form minimizers that cancel the `(i,j)` and `(j,i)` off-diagonal entries in the least-squares sense. Sweeps stop when the RMS of all `β`/`γ` of a sweep falls below `eps`, or when `itMax` sweeps have run.

---

## Requirements

- MATLAB. The code was tested on **R2015a** and uses no toolbox.
  - Exception: `simBAJD_dat` needs `mvnrnd` from the Statistics Toolbox when `condE ≠ 1`, which the tutorial uses.
- GNU Octave may work, but it has not been tested.

## Quick start

```matlab
cd CAJD
script_CAJD_tuto      % simulate data, run AJD/BAJD/CAJD, plot convergence
```

The script adds `lib/` to the path and simulates bilinear data. It then runs all three methods from the same starting point and saves the results and two figures to `figures/`:
- `CAJD_sim.mat`
- `MoMaA.png`: Moreau-Macchi index on `B`
- `MoMaE.png`: Moreau-Macchi index on `D`

> ⚠️ On macOS/Linux, change `outputDIR = '.\figures\'` to `outputDIR = fullfile('.', 'figures', filesep)` first. See [Known issues](#known-issues).

**What you should see.** Because the simulated `E` is far from orthogonal, AJD stalls far from the true solution. BAJD and CAJD both converge. With noise (SNR = 20), CAJD reaches a lower error on `B` than BAJD.

### Using it on your own data

```matlab
addpath lib
% X  : N x T x K  trials (e.g. EEG epochs, or bootstrapped ensemble averages)
% Cf : N x N x F  symmetric targets (e.g. co-spectra, covariances)

% Initialization from the SVD of the average trial (Section II-D of the paper)
[U,~,V] = svd(mean(X,3)', 0);
opt = struct('B0', inv(V)', 'D0', pinv(U)', 'eps', 1e-12, 'itMax', 100);

[B, D, S, C, info] = gp_CAJD_GPT(X, Cf, opt);

Ahat = inv(B);     % spatial patterns  (columns = scalp maps)
Ehat = pinv(D);    % temporal patterns (rows    = time courses)
% back-projection of source i for trial k:  Ahat(:,i) * S(i,i,k) * Ehat(i,:)
```

---

## Function reference (`lib/`)

| Function | Purpose |
|---|---|
| `gp_CAJD_GPT(X, C, options)` → `[B, D, S, C, info, options]` | Composite AJD of trials `X` (N×T×K) and symmetric targets `C` (N×N×F). |
| `gp_BAJD_GPT(X, options)` → `[B, D, S, info, options]` | Bilinear AJD of trials `X` only. |
| `gp_AJD_GPT(C, options)` → `[B, C, info, options]` | Classic non-orthogonal AJD of symmetric matrices `C`. |
| `simBAJD_dat(options)` → `[X, C, A, E, info]` | Simulates trials `X = A S_k Eᵀ + noise` and targets `C = A Λ_f Aᵀ + noise` that share the same `A`. |
| `convertDat_BAJD(X, Cf)` → `[C, Ct]` | Builds the AJD target set from bilinear data. `Ct(:,:,k) = cov(X_kᵀ)`, and `C = [Cf rescaled to the same total trace as Ct, Ct]`. |
| `gp_critMoAm(A, B)` | Moreau-Macchi (a.k.a. Moreau-Amari) index of `BᵀA`. It is 0 for perfect separation up to permutation and scaling. |
| `gp_critOFF(M)` | Normalized off-diagonality of a set of matrices, from 0 to 1. |
| `gp_getOpt(options, defaults)` | Merges a user `options` struct with defaults and warns on unknown fields. |
| `subplot1` | Third-party subplot helper by E. O. Ofek, used only for plotting. |

### `options` fields (GPT solvers)

| Field | Default | Meaning |
|---|---|---|
| `B0` | `eye(N)` | Initial spatial unmixing matrix (N×N). |
| `D0` | `orth(randn(T,N))` | Initial temporal unmixing matrix (T×N). Not used by AJD. |
| `eps` | `1e-12` | Stopping tolerance on the RMS update size of a sweep. |
| `itMax` | `100` | Maximum number of sweeps. |
| `A` | `[]` | True spatial mixing matrix. Only used to track `info.critA`. |
| `E` | `[]` | True temporal mixing matrix. Only used to track `info.critE`. |

### `info` (one entry per sweep, `info(1)` = initial state)

`it` (sweep number), `conv` (RMS of the sweep's β/γ), `B`, `D`, `critOFF`/`critOFFS`/`critOFFC` (off-diagonality of the transformed matrices), `critA`, `critE` (Moreau-Macchi index against the true `A`/`E`).

### `simBAJD_dat` options

| Field | Default | Meaning |
|---|---|---|
| `N`, `T` | 8, 128 | Number of sources/sensors, number of samples. |
| `K`, `F` | 20, 20 | Number of trials `X_k`, number of symmetric targets `C_f`. |
| `condA` | 1 | `1` gives an orthogonal `A`. `[lo hi]` gives a random `A` whose condition number lies in that range. |
| `condE` | 1 | Same for `E`. Use a wide range: it is found by rejection sampling. |
| `SNR` | `inf` | Noise level. Trial noise has standard deviation `1/sqrt(SNR)`. |

---

## Differences between the code and the paper

- **Weighting α.** The paper defines `f_c = α f + (1−α) f_b` and uses α = 0.5. `gp_CAJD_GPT` has no α parameter. It simply adds the bilinear and the linear sums with equal weight when computing β. Since the paper's closed form carries an extra factor 2 on the linear term, this is not exactly α = 0.5, so the balance between the two terms is set by the relative scale of `X` and `C`.
- **Tutorial vs. paper simulation.** The tutorial runs **one** small simulation: N = 8, K = 100, F = 20, orthogonal `A`, cond(`E`) in [50, 100], SNR ∈ {∞, 20}. The paper's Fig. 2 and Table I average 100 draws with N = 16 and non-orthogonal `A` and `E`. To get closer, raise `nTest`, `N` and `condA`.
- **No EEG data.** The P300 dataset (Brain Invaders) behind Figs. 3–4 is not included. The images in `pic/` are the paper's figures:
  - `EEG_independent_sources.png`: the 16 CAJD sources.
  - `evoked_sources_comparison.png`: back-projected P3b for CAJD, BAJD and AJD at AFz, Cz and Pz.

![Back-projected P3b source at AFz, Cz and Pz: raw data vs. CAJD, BAJD and AJD.](pic/evoked_sources_comparison.png)

## Known issues

These come from reading the code. They have not been re-run.

- **Windows-only output path:** `script_CAJD_tuto.m` uses `'.\figures\'`, which breaks on macOS/Linux.
- **`A`/`E` are not actually optional:** `gp_critMoAm` errors when `A` or `B` is empty because `N` is never set. The `return` on that branch is commented out. Until that is fixed, pass the true `A`/`E`, or patch `gp_critMoAm` to `return` when either input is empty.
- **`gp_CAJD_GPT(X, C)` without `options`:** this fails because the default-options check tests `nargin<2` instead of `nargin<3`. Always pass an `options` struct, even an empty `struct()`.
- **`info(1)` off-diagonality:** in `gp_AJD_GPT` and `gp_CAJD_GPT`, `info(1).critOFF*` is computed on the raw `C` rather than on `B0ᵀ C B0`.
- **Field names in BAJD:** `gp_BAJD_GPT` stores the initial off-diagonality as `critOFFS` but later sweeps as `critOFF`.

## Authors and license

Florent Bouchard and Louis Korczowski, GIPSA-lab (Univ. Grenoble Alpes, CNRS), 2016, with Christian Jutten and Marco Congedo. Supported by ERC CHESS (2012-ERC-AdG-320684) and LabEx PERSYVAL-Lab (ANR-11-LABX-0025-01).

Released under the **GNU GPL v3**, see [LICENSE](LICENSE). `lib/subplot1.m` is by Eran O. Ofek.

### Citation

```bibtex
@inproceedings{korczowski2016cajd,
  title     = {Mining the Bilinear Structure of Data with Approximate Joint Diagonalization},
  author    = {Korczowski, Louis and Bouchard, Florent and Jutten, Christian and Congedo, Marco},
  booktitle = {24th European Signal Processing Conference (EUSIPCO 2016)},
  pages     = {667--671},
  year      = {2016},
  address   = {Budapest, Hungary},
  url       = {https://hal.science/hal-01357245}
}
```
