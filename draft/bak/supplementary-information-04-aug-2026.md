# Contents of this file {.unnumbered #sec:contents}

1. Text S1
2. Figures -@fig:point-by-point-summary-est to -@fig:swp27-swp28-swp29-swp30-swp31-swp32-composition
3. Tables -@tbl:global-obs-summary-table to -@tbl:nlopt-summary

\clearpage

# Introduction {.unnumbered #sec:introduction}

This Supplementary Information (SI) file provides additional methodological details, results, and regional visualizations accompanying the main text. The SI is organized as follows.

Text S1 documents the Kriging variogram optimization procedure in full, including the ordinary Kriging equations, parameter search space, cost function, and convergence outcomes for all 34 SubMap transect-set domains.

Figures -@fig:point-by-point-summary-est to -@fig:point-by-point-summary-krg-residuals show point-by-point summary distributions of interpolation estimate differences, interpolation uncertainty differences, and interpolation residuals that support specific claims in the Results but were moved from the main text.

Figures -@fig:npa-set1-variogram to -@fig:swp-set4-variogram show the optimized experimental variogram and fitted variogram model for each of the 34 SubMap transect-set domains. These figures directly support the Kriging Optimization results in the main text, which establishes that no single variogram model or characteristic correlation length describes subduction zone surface heat flow continuity globally.

Figures -@fig:npa-set1-comp to -@fig:swp-set4-comp visualize surface heat flow data for all 34 SubMap transect sets, documenting the full range of spatial surface heat flow patterns across all analyzed domains. Figures -@fig:npa01-npa02-npa03-npa04-composition to -@fig:swp27-swp28-swp29-swp30-swp31-swp32-composition quantify the along-strike continuity of these surface heat flow patterns through CCF analysis.

Tables -@tbl:global-obs-summary-table to -@tbl:nlopt-summary contain summary statistics and cross-correlation results for all 34 transect sets, including IHFC 2024 observation statistics, Similarity--Krige estimate and uncertainty differences, interpolation residuals, global SubMap transect cross-correlations, per-transect method comparisons, and optimized Kriging parameters.

\clearpage

# Definition of Symbols {.unnumbered #sec:symbols}

|Parameter                                      |Symbol                        |Unit             |Equations                                          |
|:----------------------------------------------|:-----------------------------|:----------------|:--------------------------------------------------|
|**Spatial Field**                              |                              |                 |                                                   |
|Spatial location                               |$u$                           |km               |1, [-@eq:krg-assump; -@eq:blue; -@eq:min-var]      |
|Observation (surface heat flow) at $u$         |$Z(u)$                        |mW m$^{-2}$      |1, [-@eq:krg-assump; -@eq:min-var]                 |
|Mean surface heat flow                         |$\mu$                         |mW m$^{-2}$      |[-@eq:krg-assump]                                  |
|Covariance function                            |$C(h)$                        |mW$^2$ m$^{-4}$  |[-@eq:krg-assump; -@eq:min-var]                    |
|Raw sample count                               |$N_\mathrm{obs}$              |-                |-                                                  |
|**Variogram Model**                            |                              |                 |                                                   |
|Separation distance (lag) between observations |$h$                           |km               |1, 2, [-@eq:krg-assump; -@eq:vgrm-mods, -@eq:bins] |
|Maximum separation distance                    |$\max(h)$                     |km               |2, [-@eq:bins]                                     |
|Lag cutoff constant                            |$\kappa$                      |-                |2, [-@eq:bins]                                     |
|No. of lag intervals                           |$n$                           |-                |2, [-@eq:bins; -@eq:cost]                          |
|Lag bin width                                  |$\delta$                      |km               |[-@eq:bins]                                        |
|Observation pair count at lag $h$              |$N(h)$                        |-                |1, 2                                               |
|Experimental variogram                         |$\hat{\gamma}(h)$             |mW$^2$ m$^{-4}$  |1                                                  |
|Variogram model                                |$\gamma(h)$                   |mW$^2$ m$^{-4}$  |1, [-@eq:vgrm-mods]                                |
|Variogram model type                           |Sph, Exp, Gau                 |-                |[-@eq:vgrm-mods]                                   |
|Nugget variance                                |$\tau^2$                      |mW$^2$ m$^{-4}$  |[-@eq:vgrm-mods]                                   |
|Sill variance                                  |$\sigma^2$                    |mW$^2$ m$^{-4}$  |[-@eq:vgrm-mods]                                   |
|Effective range                                |$a$                           |km               |[-@eq:vgrm-mods]                                   |
|**Kriging Estimation**                         |                              |                 |                                                   |
|Krige estimate at $u$                          |$\hat{Z}(u)$                  |mW m$^{-2}$      |[-@eq:blue; -@eq:min-var]                          |
|No. of observations used for Kriging           |$n_K$                         |-                |[-@eq:blue; -@eq:min-var]                          |
|Kriging weight                                 |$\lambda_i$                   |-                |[-@eq:blue; -@eq:min-var]                          |
|Kriging variance                               |$\sigma_K^2(u)$               |mW$^2$ m$^{-4}$  |[-@eq:min-var]                                     |
|Maximum neighborhood observations              |$m$                           |-                |-                                                  |
|Maximum neighborhood size                      |$d_\mathrm{max}$              |km               |-                                                  |
|**Optimization**                               |                              |                 |                                                   |
|Optimization objective                         |$J$                           |-                |[-@eq:cost]                                        |
|Variogram weighted least squares error         |$\mathrm{WLS}_\mathrm{error}$ | mW$^4$ m$^{-8}$ |[-@eq:cost]                                        |
|Standard deviation of experimental variogram   |$\sigma_{\gamma}$             | mW$^2$ m$^{-4}$ |[-@eq:cost]                                        |
|Kriging Root Mean Square Error                 |$\mathrm{RMSE}_\mathrm{Krg}$  | mW m$^{-2}$     |[-@eq:cost]                                        |
|Standard deviation of cross-validation         |$\sigma_\mathrm{cv}$          | mW m$^{-2}$     |[-@eq:cost]                                        |

\clearpage

# Text S1. Kriging Variogram Optimization {.unnumbered #sec:kriging-system-and-optimization}

## Overview {.unnumbered #sec:overview}

Kriging was applied to IHFC 2024 heat flow observations within each of 34 SubMap transect-set domains. For each domain, the variogram model type and all associated parameters were selected by constrained nonlinear optimization rather than by manual fitting, following the approach of @li2018. The optimization minimized a composite cost function combining variogram-fitting error and cross-validation interpolation error. This section documents the optimization procedure, parameter search space, cost function, and convergence results for all 34 domains.

## Ordinary Kriging Equations {.unnumbered #sec:ordinary-kriging-equations}

This study applies local isotropic ordinary Kriging methods under the following general assumptions:

- $\hat{\gamma}(h)$ is directionally invariant (isotropic)
- $\hat{\gamma}(h)$ is evaluated in two-dimensions and neglects elevation
- The first and second moments of $Z(u)$ are assumed to follow the conditions:

$$
  \begin{aligned}
    &E[Z(u)] = \mu = constant \\
    &E[(Z(u + h) - \mu)(Z(u) - \mu)] = C(h)
  \end{aligned}
$$ {#eq:krg-assump}

where $h$ is the lag distance, $C(h)$ is the covariance function, $E[Z(u)]$ is the expected value of the random variable $Z(u)$, and $\mu$ is the arithmetic mean of $Z(u)$.

Equation @eq:krg-assump is known as "weak second-order stationarity". It assumes the underlying probability distribution of the observations $Z(u)$ does not change in space and the covariance $C(h)$ only depends on the distance $h$ between two observations. These assumptions are expected to be valid in cases where the underlying natural process is stochastic, spatially continuous, and has the property of additivity such that $\frac{1}{n_K}\sum_{i=1}^{n_K} Z(u_i)$ has the same meaning as $Z(u)$ [@bardossy1997].

Equation @eq:krg-assump is likely invalid in regions that transition between two or more tectonic regimes, however. For example, the expected (mean) surface heat flow $E[Z(u)]$ will change when moving from a spreading center to a subduction zone and thus $E[Z(u)] \neq constant$ over the region of interest. In other words, stationarity is violated and Kriging estimates may become spurious. Careful selection of Kriging parameters (outlined below; e.g. maximum point-pairs to use for local Kriging) can reduce or eliminate violations of stationarity assumptions embodied in @eq:krg-assump.

The second step is fitting a variogram model $\gamma(h)$ to the experimental variogram $\hat{\gamma}(h)$. This study fits three popular variogram models with sills $\sigma^2$ to the experimental variogram. The models are defined as [@pebesma2004]:

$$
  \begin{aligned}
    Sph &\leftarrow \gamma(h) =
    \begin{cases}
      \tau^2 + \sigma^2 \left[ \frac{3}{2}\, \frac{h}{a} - \frac{1}{2}\, \left( \frac{h}{a} \right)^3 \right] \quad \text{for } \  0 \leq h \leq a \\
      \tau^2 + \sigma^2 \quad \text{for } \  h > a
    \end{cases} \\
    Exp &\leftarrow \gamma(h) = \tau^2 + \sigma^2 \left[ 1 - \exp(\frac{-h}{a}) \right] \quad \text{for } \  h \geq 0 \\
    Gau &\leftarrow \gamma(h) = \tau^2 + \sigma^2 \left[ 1 - \exp(\left(\frac{-h}{a}\right)^2) \right] \quad \text{for } \  h \geq 0 \\
  \end{aligned}
$$ {#eq:vgrm-mods}

where $h$ is the lag distance, $\tau^2$ is the nugget, $\sigma^2$ is the sill, $a$ is the effective range. The models are Spherical, Exponential, and Gaussian. For models without explicit sills (Exp and Gau), the effective range $a$ is the distance where the variogram reaches 95% of its maximum defined as 3$a$ and $\sqrt{3}a$ for Exp and Gau, respectively [@graler2016; @pebesma2004]. The function `fit.variogram` in `gstat` is used to try all variogram models. The best model is selected by the minimum weighted least squares error [@pebesma2004] with weights proportional to the number of points in each lag divided by the squared lag distance $N(h_k)\, /\, h_k^2$.

Ordinary Kriging is used for interpolation, which estimates unknown observations $\hat{Z}(u)$ as a linear combination of all known observations [@bardossy1997]:

$$
  \hat{Z}(u) = \sum_{i=1}^{n_K} \lambda_i Z(u_i)
$$ {#eq:blue}

The conditions in Equation @eq:krg-assump set up a constrained minimization problem that can be solved with a system of linear equations. The expected value of $Z(u)$ is assumed to be the mean according to @eq:krg-assump, so the weights must be:

$$
  \begin{aligned}
    E[\hat{Z}(u)] &= \sum_{i=1}^{n_K} \lambda_i E[Z(u_i)] \\
    \sum_{i=1}^{n_K} \lambda_i &= 1
  \end{aligned}
$$ {#eq:unbiased}

This constraint is known as the unbiased condition, which states that the sum of the weights must equal one. However, there is an infinite set of real numbers one could use for the weights, $\lambda_i$. The goal is to find the set of weights in Equation @eq:blue that minimizes the estimation variance. This can be solved by minimizing the covariance function, $C(h)$ from Equation @eq:krg-assump:

$$
  \begin{aligned}
    & \sigma_K^2(u) = Var[Z(u) - \hat{Z}(u)] = \\
    & E\left[(Z(u) - \sum_{i=1}^{n_K} \lambda_i Z(u_i))^2\right] = \\
    & E\left[Z(u)^2 + \sum_{j=1}^{n_K} \sum_{i=1}^{n_K} \lambda_j \lambda_i Z(u_j)Z(u_i) - 2 \sum_{i=1}^{n_K} \lambda_i Z(u_i)Z(u)\right] = \\
    & C(0) + \sum_{j=1}^{n_K} \sum_{i=1}^{n_K} \lambda_j \lambda_i C(u_i - u_j) - 2 \sum_{i=1}^{n_K} \lambda_i C(u_i - u)
  \end{aligned}
$$ {#eq:min-var}

Minimizing Equation @eq:min-var subject to the unbiased condition (Equation @eq:unbiased), yields the best linear unbiased estimator [BLUE, @bardossy1997] for Equation @eq:blue. Together Equations @eq:blue and @eq:unbiased comprise the ordinary Kriging system. The functions `krige` and `krige.cv` in `gstat` are used for surface heat flow interpolation and error estimation by k-fold cross-validation [@pebesma2004].

## Parameter Search Space {.unnumbered #sec:parameter-search-space}

Four parameters were optimized for each domain:

1. *Number of lag intervals ($n$)*: integer values on [20, 50]
2. *Variogram cutoff constant ($\kappa$):* values on [3, 12]; limits the maximum lag used to construct the experimental variogram to $\max(h)\, /\, \kappa$
3. *Maximum neighborhood observations ($m$):* integer values on [2, 50]; the number of nearest observations used for each estimate
4. *Maximum neighborhood size ($d_\mathrm{max}$):* values on [100, 500] km; the maximum distance for including observations in Kriging

The bin width $\delta$ was computed as:

$$
  \delta = \frac{\max(h)}{\kappa\, n}
$$ {#eq:bins}

where $\max(h)$ is the maximum pairwise separation distance within the convex polygon domain enclosing transect sets.

## Cost Function {.unnumbered #sec:cost-function}

The composite cost was defined as:

$$
  J = \frac{(\mathrm{WLS}_\mathrm{error}\, /\, n)^{1/4}}{2\, \sigma_{\gamma}^{1/2}}  + \frac{\mathrm{RMSE}_\mathrm{Krg}}{2\, \sigma_\mathrm{cv}}
$$ {#eq:cost}

where $J$ is the optimization objective to minimize, $(\mathrm{WLS}_\mathrm{error}\, /\, n)^{1/4}$ quantified the weighted deviation of the fitted parametric variogram model from the experimental variogram $\hat{\gamma}(h)$ [@pebesma2004], and $\mathrm{RMSE}_\mathrm{Krg}$ quantified leave-one-out cross-validation error computed over the training observations within the domain. Both components were normalized by their respective standard deviations to produce dimensionless, comparably scaled terms. Equal weighting (0.5) was applied to both terms.

## Optimization Algorithm {.unnumbered #sec:optimization-algorithm}

To ensure a robust exploration of the parameter space and mitigate the risk of local minima, a two-stage global-to-local optimization workflow was implemented via the R package `nloptr`.

First, an initial exploratory phase was conducted across the four-dimensional continuous parameter space using Latin Hypercube Sampling (LHS). For each model type, 100 starting configurations were chosen by LHS across the predefined parameter bounds using a fixed random seed (42). The cost function was evaluated at each LHS candidate point, and the highest-performing parameter configurations (those yielding the lowest total costs) were selected to serve as initial coordinates for the formal optimization.

Global optimization was then initiated from these best-performing LHS candidates using the Multi-Level Single-Linkage (MLSL) algorithm, paired with a Nelder-Mead local solver. Following the global MLSL search, a second fine-tuning pass was performed, which searched a restricted neighborhood around the best solution found in the first pass. The model type (Exp, Sph, Gau) was treated as a discrete outer parameter: the continuous LHS-driven optimization sequence was run independently for each model type, and the model yielding the lowest total cost was selected.

## Convergence Summary {.unnumbered #sec:convergence-summary}

All 34 optimization runs converged successfully. For each model type within a domain, the parameter space was initially evaluated using 100 LHS exploratory samples, from which the 3 best-performing starts were selected for refinement. Formal optimization converged to a solution within 100 MLSL evaluations per start, followed by up to 500 Nelder-Mead fine-tuning evaluations.

Convergence was defined as a change in total cost below 1e-08 between successive iterations in the fine-tuning pass. The selected model type, converged parameter values, and total cost for each domain are summarized in Table @tbl:nlopt-summary. Variogram plots showing the experimental variogram points and fitted model curve for each domain are provided in Figures -@fig:npa-set1-variogram to -@fig:swp-set4-variogram.

No systematic relationship was observed between the number of optimization iterations and the total cost, or between observation count and the selected model type. Domains with comparable observation counts frequently converged on different model types and spatial scales, and similarly low-cost solutions arose from a wide range of parameter configurations. These outcomes confirm that spatial continuity in subduction zone surface heat flow cannot be characterized by a single variogram form or correlation length.

\clearpage

# Supplementary Maps and Surface Heat Flow Profiles (Figures S1--S4) {.unnumbered #sec:supplementary-maps-and-heat-flow-profiles}

|Label     |Description                            |
|:---------|:--------------------------------------|
|Figure S1 |Point-by-point estimate differences    |
|Figure S2 |Point-by-point uncertainty differences |
|Figure S3 |Point-by-point Sim--Obs residuals      |
|Figure S4 |Point-by-point Krg--Obs residuals      |

\clearpage

![Point-by-point distributions of Similarity--Krige surface heat flow estimate differences (mW m$^{-2}$) within each SubMap transect-set domain. Distributions are colored by quartiles (25%, 50%, 75%). NPA: North Pacific--Aleutians; SAM: South America; SEA: Southeast Asia; SWP: Southwest Pacific. See @tbl:est-diff-summary for summary statistics.](../figs/summary/point-by-point-summary-est.png){#fig:point-by-point-summary-est}

![Point-by-point distributions of Similarity--Krige estimate-uncertainty differences (mW m$^{-2}$) within each SubMap transect-set domain. Distributions are colored by quartiles (25%, 50%, 75%). NPA: North Pacific--Aleutians; SAM: South America; SEA: Southeast Asia; SWP: Southwest Pacific. See @tbl:sigma-diff-summary for summary statistics.](../figs/summary/point-by-point-summary-sigma.png){#fig:point-by-point-summary-sigma}

![Point-by-point distributions of Sim--Obs residuals (mW m$^{-2}$) for each SubMap transect set. Distributions are colored by quartiles (25%, 50%, 75%). NPA: North Pacific--Aleutians; SAM: South America; SEA: Southeast Asia; SWP: Southwest Pacific. See @tbl:sim-residual-summary for summary statistics.](../figs/summary/point-by-point-summary-sim-residuals.png){#fig:point-by-point-summary-sim-residuals}

![Point-by-point distributions of Krg--Obs residuals (mW m$^{-2}$) for each SubMap transect set. Distributions are colored by quartiles (25%, 50%, 75%). NPA: North Pacific--Aleutians; SAM: South America; SEA: Southeast Asia; SWP: Southwest Pacific. See @tbl:krg-residual-summary for summary statistics.](../figs/summary/point-by-point-summary-krg-residuals.png){#fig:point-by-point-summary-krg-residuals}

\clearpage

# Optimized Experimental Variogram (Figures S5--S38) {.unnumbered #sec:optimized-experimental-variogram}

|Label      |Description                                                |
|:----------|:----------------------------------------------------------|
|Figure S5  |Variogram and model for NPA_SET1 (Japan)                   |
|Figure S6  |Variogram and model for NPA_SET2 (Kurils)                  |
|Figure S7  |Variogram and model for NPA_SET3 (Kurils--Kamchatka)       |
|Figure S8  |Variogram and model for NPA_SET4 (Aleutians)               |
|Figure S9  |Variogram and model for NPA_SET5 (Aleutians)               |
|Figure S10 |Variogram and model for NPA_SET6 (Alaska)                  |
|Figure S11 |Variogram and model for NPA_SET7 (Alaska)                  |
|Figure S12 |Variogram and model for NPA_SET8 (Cascadia)                |
|Figure S13 |Variogram and model for SAM_SET1 (Central America)         |
|Figure S14 |Variogram and model for SAM_SET2 (Central America--Panama) |
|Figure S15 |Variogram and model for SAM_SET3 (Venezuela)               |
|Figure S16 |Variogram and model for SAM_SET4 (Muertos--Antilles)       |
|Figure S17 |Variogram and model for SAM_SET5 (Andean)                  |
|Figure S18 |Variogram and model for SAM_SET6 (Andean)                  |
|Figure S19 |Variogram and model for SAM_SET7 (Andean)                  |
|Figure S20 |Variogram and model for SAM_SET8 (Andean)                  |
|Figure S21 |Variogram and model for SAM_SET9 (Andean)                  |
|Figure S22 |Variogram and model for SAM_SET10 (Andean)                 |
|Figure S23 |Variogram and model for SEA_SET1 (Andaman--Sumatra)        |
|Figure S24 |Variogram and model for SEA_SET2 (Sumatra)                 |
|Figure S25 |Variogram and model for SEA_SET3 (Java)                    |
|Figure S26 |Variogram and model for SEA_SET4 (Java--Flores)            |
|Figure S27 |Variogram and model for SEA_SET5 (Timor--Seram)            |
|Figure S28 |Variogram and model for SEA_SET6 (Halmahera--Manila)       |
|Figure S29 |Variogram and model for SEA_SET7 (Manila)                  |
|Figure S30 |Variogram and model for SEA_SET8 (Ryukyus)                 |
|Figure S31 |Variogram and model for SEA_SET9 (Ryukyus--Sagami)         |
|Figure S32 |Variogram and model for SEA_SET10 (Izu--Bonin)             |
|Figure S33 |Variogram and model for SEA_SET11 (Mariana)                |
|Figure S34 |Variogram and model for SEA_SET12 (Yap--Palau)             |
|Figure S35 |Variogram and model for SWP_SET1 (New Britain--Solomons)   |
|Figure S36 |Variogram and model for SWP_SET2 (New Hebrides)            |
|Figure S37 |Variogram and model for SWP_SET3 (Tonga--Kermadec)         |
|Figure S38 |Variogram and model for SWP_SET4 (Kermadec--Hikurangi)     |

\clearpage

![Optimized experimental variogram (points) and fitted exponential model (curve) for NPA_SET1 (Japan). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/NPA_SET1-variogram.png){#fig:npa-set1-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for NPA_SET2 (Kurils). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/NPA_SET2-variogram.png){#fig:npa-set2-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for NPA_SET3 (Kurils--Kamchatka). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/NPA_SET3-variogram.png){#fig:npa-set3-variogram}

![Optimized experimental variogram (points) and fitted Gaussian model (curve) for NPA_SET4 (Aleutians). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/NPA_SET4-variogram.png){#fig:npa-set4-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for NPA_SET5 (Aleutians). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/NPA_SET5-variogram.png){#fig:npa-set5-variogram}

![Optimized experimental variogram (points) and fitted Gaussian model (curve) for NPA_SET6 (Alaska). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/NPA_SET6-variogram.png){#fig:npa-set6-variogram}

![Optimized experimental variogram (points) and fitted spherical model (curve) for NPA_SET7 (Alaska). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/NPA_SET7-variogram.png){#fig:npa-set7-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for NPA_SET8 (Cascadia). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/NPA_SET8-variogram.png){#fig:npa-set8-variogram}

![Optimized experimental variogram (points) and fitted Gaussian model (curve) for SAM_SET1 (Central America). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SAM_SET1-variogram.png){#fig:sam-set1-variogram}

![Optimized experimental variogram (points) and fitted spherical model (curve) for SAM_SET2 (Central America--Panama). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SAM_SET2-variogram.png){#fig:sam-set2-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for SAM_SET3 (Venezuela). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SAM_SET3-variogram.png){#fig:sam-set3-variogram}

![Optimized experimental variogram (points) and fitted Gaussian model (curve) for SAM_SET4 (Muertos--Antilles). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SAM_SET4-variogram.png){#fig:sam-set4-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for SAM_SET5 (Andean). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SAM_SET5-variogram.png){#fig:sam-set5-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for SAM_SET6 (Andean). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SAM_SET6-variogram.png){#fig:sam-set6-variogram}

![Optimized experimental variogram (points) and fitted spherical model (curve) for SAM_SET7 (Andean). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SAM_SET7-variogram.png){#fig:sam-set7-variogram}

![Optimized experimental variogram (points) and fitted Gaussian model (curve) for SAM_SET8 (Andean). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SAM_SET8-variogram.png){#fig:sam-set8-variogram}

![Optimized experimental variogram (points) and fitted spherical model (curve) for SAM_SET9 (Andean). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SAM_SET9-variogram.png){#fig:sam-set9-variogram}

![Optimized experimental variogram (points) and fitted spherical model (curve) for SAM_SET10 (Andean). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SAM_SET10-variogram.png){#fig:sam-set10-variogram}

![Optimized experimental variogram (points) and fitted spherical model (curve) for SEA_SET1 (Andaman--Sumatra). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SEA_SET1-variogram.png){#fig:sea-set1-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for SEA_SET2 (Sumatra). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SEA_SET2-variogram.png){#fig:sea-set2-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for SEA_SET3 (Java). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SEA_SET3-variogram.png){#fig:sea-set3-variogram}

![Optimized experimental variogram (points) and fitted spherical model (curve) for SEA_SET4 (Java--Flores). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SEA_SET4-variogram.png){#fig:sea-set4-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for SEA_SET5 (Timor--Seram). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SEA_SET5-variogram.png){#fig:sea-set5-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for SEA_SET6 (Halmahera--Manila). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SEA_SET6-variogram.png){#fig:sea-set6-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for SEA_SET7 (Manila). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SEA_SET7-variogram.png){#fig:sea-set7-variogram}

![Optimized experimental variogram (points) and fitted spherical model (curve) for SEA_SET8 (Ryukyus). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SEA_SET8-variogram.png){#fig:sea-set8-variogram}

![Optimized experimental variogram (points) and fitted spherical model (curve) for SEA_SET9 (Ryukyus--Sagami). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SEA_SET9-variogram.png){#fig:sea-set9-variogram}

![Optimized experimental variogram (points) and fitted spherical model (curve) for SEA_SET10 (Izu--Bonin). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SEA_SET10-variogram.png){#fig:sea-set10-variogram}

![Optimized experimental variogram (points) and fitted Gaussian model (curve) for SEA_SET11 (Mariana). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SEA_SET11-variogram.png){#fig:sea-set11-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for SEA_SET12 (Yap--Palau). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SEA_SET12-variogram.png){#fig:sea-set12-variogram}

![Optimized experimental variogram (points) and fitted exponential model (curve) for SWP_SET1 (New Britain--Solomons). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SWP_SET1-variogram.png){#fig:swp-set1-variogram}

![Optimized experimental variogram (points) and fitted spherical model (curve) for SWP_SET2 (New Hebrides). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SWP_SET2-variogram.png){#fig:swp-set2-variogram}

![Optimized experimental variogram (points) and fitted Gaussian model (curve) for SWP_SET3 (Tonga--Kermadec). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SWP_SET3-variogram.png){#fig:swp-set3-variogram}

![Optimized experimental variogram (points) and fitted Gaussian model (curve) for SWP_SET4 (Kermadec--Hikurangi). Variogram parameters and optimization cost are given in Table @tbl:nlopt-summary.](../figs/submap/optimal/SWP_SET4-variogram.png){#fig:swp-set4-variogram}

\clearpage

# Regional Surface Heat Flow Comparisons (Figures S39--S71) {.unnumbered #sec:regional-heat-flow-comparisons}

|Label      |Description                                               |
|:----------|:---------------------------------------------------------|
|Figure S39 |Spatial comparison for NPA_SET1 (Japan)                   |
|Figure S40 |Spatial comparison for NPA_SET3 (Kurils--Kamchatka)       |
|Figure S41 |Spatial comparison for NPA_SET4 (Aleutians)               |
|Figure S42 |Spatial comparison for NPA_SET5 (Aleutians)               |
|Figure S43 |Spatial comparison for NPA_SET6 (Alaska)                  |
|Figure S44 |Spatial comparison for NPA_SET7 (Alaska)                  |
|Figure S45 |Spatial comparison for NPA_SET8 (Cascadia)                |
|Figure S46 |Spatial comparison for SAM_SET2 (Central America)         |
|Figure S47 |Spatial comparison for SAM_SET2 (Central America--Panama) |
|Figure S48 |Spatial comparison for SAM_SET3 (Venezuela)               |
|Figure S49 |Spatial comparison for SAM_SET4 (Muertos--Antilles)       |
|Figure S50 |Spatial comparison for SAM_SET5 (Andean)                  |
|Figure S51 |Spatial comparison for SAM_SET6 (Andean)                  |
|Figure S52 |Spatial comparison for SAM_SET7 (Andean)                  |
|Figure S53 |Spatial comparison for SAM_SET8 (Andean)                  |
|Figure S54 |Spatial comparison for SAM_SET9 (Andean)                  |
|Figure S55 |Spatial comparison for SAM_SET10 (Andean)                 |
|Figure S56 |Spatial comparison for SEA_SET1 (Andaman--Sumatra)        |
|Figure S57 |Spatial comparison for SEA_SET2 (Sumatra)                 |
|Figure S58 |Spatial comparison for SEA_SET2 (Java)                    |
|Figure S59 |Spatial comparison for SEA_SET4 (Java--Flores)            |
|Figure S60 |Spatial comparison for SEA_SET5 (Timor--Seram)            |
|Figure S61 |Spatial comparison for SEA_SET6 (Halmahera--Manila)       |
|Figure S62 |Spatial comparison for SEA_SET7 (Manila)                  |
|Figure S63 |Spatial comparison for SEA_SET8 (Ryukyus)                 |
|Figure S64 |Spatial comparison for SEA_SET9 (Ryukyus--Sagami)         |
|Figure S65 |Spatial comparison for SEA_SET10 (Izu--Bonin)             |
|Figure S66 |Spatial comparison for SEA_SET11 (Mariana)                |
|Figure S67 |Spatial comparison for SEA_SET12 (Yap--Palau)             |
|Figure S68 |Spatial comparison for SWP_SET1 (New Britain--Solomons)   |
|Figure S69 |Spatial comparison for SWP_SET2 (New Hebrides)            |
|Figure S70 |Spatial comparison for SWP_SET3 (Tonga--Kermadec)         |
|Figure S71 |Spatial comparison for SWP_SET4 (Kermadec--Hikurangi)     |

\clearpage

![Surface heat flow observations, Krige estimates, and Similarity estimates for NPA_SET1 (Japan). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/NPA_SET1-comp.png){#fig:npa-set1-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for NPA_SET3 (Kurils--Kamchatka). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/NPA_SET3-comp.png){#fig:npa-set3-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for NPA_SET4 (Aleutians). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/NPA_SET4-comp.png){#fig:npa-set4-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for NPA_SET4 (Aleutians). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/NPA_SET5-comp.png){#fig:npa-set5-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for NPA_SET6 (Alaska). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/NPA_SET6-comp.png){#fig:npa-set6-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for NPA_SET7 (Alaska). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/NPA_SET7-comp.png){#fig:npa-set7-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for NPA_SET7 (Cascadia). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/NPA_SET8-comp.png){#fig:npa-set8-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SAM_SET2 (Central America). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SAM_SET1-comp.png){#fig:sam-set1-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SAM_SET2 (Central America--Panama). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SAM_SET2-comp.png){#fig:sam-set2-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SAM_SET3 (Venezuela). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SAM_SET3-comp.png){#fig:sam-set3-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SAM_SET4 (Muertos--Antilles). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SAM_SET4-comp.png){#fig:sam-set4-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SAM_SET5 (Andean). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SAM_SET5-comp.png){#fig:sam-set5-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SAM_SET6 (Andean). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SAM_SET6-comp.png){#fig:sam-set6-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SAM_SET7 (Andean). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SAM_SET7-comp.png){#fig:sam-set7-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SAM_SET8 (Andean). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SAM_SET8-comp.png){#fig:sam-set8-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SAM_SET9 (Andean). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SAM_SET9-comp.png){#fig:sam-set9-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SAM_SET10 (Andean). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SAM_SET10-comp.png){#fig:sam-set10-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SEA_SET1 (Andaman--Sumatra). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SEA_SET1-comp.png){#fig:sea-set1-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SEA_SET2 (Sumatra). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SEA_SET2-comp.png){#fig:sea-set2-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SEA_SET4 (Java). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SEA_SET3-comp.png){#fig:sea-set3-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SEA_SET4 (Java--Flores). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SEA_SET4-comp.png){#fig:sea-set4-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SEA_SET5 (Timor--Seram). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SEA_SET5-comp.png){#fig:sea-set5-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SEA_SET6 (Halmahera--Manila). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SEA_SET6-comp.png){#fig:sea-set6-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SEA_SET7 (Manila). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SEA_SET7-comp.png){#fig:sea-set7-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SEA_SET8 (Ryukyus). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SEA_SET8-comp.png){#fig:sea-set8-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SEA_SET9 (Ryukyus--Sagami). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SEA_SET9-comp.png){#fig:sea-set9-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SEA_SET10 (Izu--Bonin). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SEA_SET10-comp.png){#fig:sea-set10-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SEA_SET11 (Mariana). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SEA_SET11-comp.png){#fig:sea-set11-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SEA_SET12 (Yap--Palau). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SEA_SET12-comp.png){#fig:sea-set12-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SWP_SET1 (New Britain--Solomons). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SWP_SET1-comp.png){#fig:swp-set1-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SWP_SET2 (New Hebrides). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SWP_SET2-comp.png){#fig:swp-set2-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SWP_SET3 (Tonga--Kermadec). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SWP_SET3-comp.png){#fig:swp-set3-comp}

![Surface heat flow observations, Krige estimates, and Similarity estimates for SWP_SET4 (Kermadec--Hikurangi). Annotations are as described in Figure 2 of the main text.](../figs/submap/composition/SWP_SET4-comp.png){#fig:swp-set4-comp}

\clearpage

# CCF Analysis of Neighboring Transects (Figures S72--S122) {.unnumbered #sec:ccf-analysis-of-neighboring-transects}

|Label       |Description                                |
|:-----------|:------------------------------------------|
|Figure S72  |CCF analysis and profiles for NPA01--NPA04 |
|Figure S73  |CCF analysis and profiles for NPA05--NPA10 |
|Figure S74  |CCF analysis and profiles for NPA11--NPA14 |
|Figure S75  |CCF analysis and profiles for NPA15--NPA19 |
|Figure S76  |CCF analysis and profiles for NPA20--NPA22 |
|Figure S77  |CCF analysis and profiles for NPA23--NPA27 |
|Figure S78  |CCF analysis and profiles for NPA28--NPA30 |
|Figure S79  |CCF analysis and profiles for NPA31--NPA36 |
|Figure S80  |CCF analysis and profiles for NPA37--NPA42 |
|Figure S81  |CCF analysis and profiles for NPA43--NPA47 |
|Figure S82  |CCF analysis and profiles for SAM01--SAM05 |
|Figure S83  |CCF analysis and profiles for SAM06--SAM08 |
|Figure S84  |CCF analysis and profiles for SAM09--SAM13 |
|Figure S85  |CCF analysis and profiles for SAM14--SAM15 |
|Figure S86  |CCF analysis and profiles for SAM16--SAM21 |
|Figure S87  |CCF analysis and profiles for SAM22--SAM26 |
|Figure S88  |CCF analysis and profiles for SAM27--SAM31 |
|Figure S89  |CCF analysis and profiles for SAM32--SAM33 |
|Figure S90  |CCF analysis and profiles for SAM34--SAM38 |
|Figure S91  |CCF analysis and profiles for SAM39--SAM40 |
|Figure S92  |CCF analysis and profiles for SAM41--SAM46 |
|Figure S93  |CCF analysis and profiles for SAM47--SAM51 |
|Figure S94  |CCF analysis and profiles for SAM52--SAM56 |
|Figure S95  |CCF analysis and profiles for SAM57--SAM61 |
|Figure S96  |CCF analysis and profiles for SAM62--SAM66 |
|Figure S97  |CCF analysis and profiles for SEA01--SEA05 |
|Figure S98  |CCF analysis and profiles for SEA06--SEA07 |
|Figure S99  |CCF analysis and profiles for SEA08--SEA12 |
|Figure S100 |CCF analysis and profiles for SEA13--SEA17 |
|Figure S101 |CCF analysis and profiles for SEA18--SEA22 |
|Figure S102 |CCF analysis and profiles for SEA23--SEA27 |
|Figure S103 |CCF analysis and profiles for SEA28--SEA34 |
|Figure S104 |CCF analysis and profiles for SEA29--SEA30 |
|Figure S105 |CCF analysis and profiles for SEA35--SEA39 |
|Figure S106 |CCF analysis and profiles for SEA40--SEA44 |
|Figure S107 |CCF analysis and profiles for SEA45--SEA49 |
|Figure S108 |CCF analysis and profiles for SEA50--SEA54 |
|Figure S109 |CCF analysis and profiles for SEA55--SEA60 |
|Figure S110 |CCF analysis and profiles for SEA61--SEA66 |
|Figure S111 |CCF analysis and profiles for SEA67--SEA72 |
|Figure S112 |CCF analysis and profiles for SEA73--SEA77 |
|Figure S113 |CCF analysis and profiles for SEA78--SEA82 |
|Figure S114 |CCF analysis and profiles for SEA83--SEA86 |
|Figure S115 |CCF analysis and profiles for SEA87--SEA90 |
|Figure S116 |CCF analysis and profiles for SWP01--SWP05 |
|Figure S117 |CCF analysis and profiles for SWP06--SWP09 |
|Figure S118 |CCF analysis and profiles for SWP10--SWP14 |
|Figure S119 |CCF analysis and profiles for SWP15--SWP16 |
|Figure S120 |CCF analysis and profiles for SWP17--SWP21 |
|Figure S121 |CCF analysis and profiles for SWP22--SWP26 |
|Figure S122 |CCF analysis and profiles for SWP27--SWP32 |

\clearpage

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects NPA01--NPA04.](../figs/submap/case_study/NPA01-NPA02-NPA03-NPA04-composition.png){#fig:npa01-npa02-npa03-npa04-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects NPA05--NPA10.](../figs/submap/case_study/NPA05-NPA06-NPA07-NPA08-NPA09-NPA10-composition.png){#fig:npa05-npa06-npa07-npa08-npa09-npa10-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects NPA11--NPA14.](../figs/submap/case_study/NPA11-NPA12-NPA13-NPA14-composition.png){#fig:npa11-npa12-npa13-npa14-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects NPA15--NPA19.](../figs/submap/case_study/NPA15-NPA16-NPA17-NPA18-NPA19-composition-part1.png){#fig:npa15-npa16-npa17-npa18-npa19-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects NPA20--NPA22.](../figs/submap/case_study/NPA20-NPA21-NPA22-composition-part2.png){#fig:npa20-npa21-npa22-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring Aleutian transects NPA23--NPA27.](../figs/submap/case_study/NPA23-NPA24-NPA25-NPA26-NPA27-composition-part1.png){#fig:npa23-npa24-npa25-npa26-npa27-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects NPA28--NPA30.](../figs/submap/case_study/NPA28-NPA29-NPA30-composition-part2.png){#fig:npa28-npa29-npa30-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects NPA31--NPA36.](../figs/submap/case_study/NPA31-NPA32-NPA33-NPA34-NPA35-NPA36-composition.png){#fig:npa31-npa32-npa33-npa34-npa35-npa36-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects NPA37--NPA42.](../figs/submap/case_study/NPA37-NPA38-NPA39-NPA40-NPA41-NPA42-composition.png){#fig:npa37-npa38-npa39-npa40-npa41-npa42-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects NPA43--NPA47.](../figs/submap/case_study/NPA43-NPA44-NPA45-NPA46-NPA47-composition.png){#fig:npa43-npa44-npa45-npa46-npa47-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM01--SAM05.](../figs/submap/case_study/SAM01-SAM02-SAM03-SAM04-SAM05-composition-part1.png){#fig:sam01-sam02-sam03-sam04-sam05-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM06--SAM08.](../figs/submap/case_study/SAM06-SAM07-SAM08-composition-part2.png){#fig:sam06-sam07-sam08-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM09--SAM13.](../figs/submap/case_study/SAM09-SAM10-SAM11-SAM12-SAM13-composition-part1.png){#fig:sam09-sam10-sam11-sam12-sam13-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM14--SAM15.](../figs/submap/case_study/SAM14-SAM15-composition-part2.png){#fig:sam14-sam15-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM16--SAM21.](../figs/submap/case_study/SAM16-SAM17-SAM18-SAM19-SAM20-SAM21-composition.png){#fig:sam16-sam17-sam18-sam19-sam20-sam21-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM22--SAM26.](../figs/submap/case_study/SAM22-SAM23-SAM24-SAM25-SAM26-composition-part1.png){#fig:sam22-sam23-sam24-sam25-sam26-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM27--SAM31.](../figs/submap/case_study/SAM27-SAM28-SAM29-SAM30-SAM31-composition-part2.png){#fig:sam27-sam28-sam29-sam30-sam31-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM32--SAM33.](../figs/submap/case_study/SAM32-SAM33-composition-part3.png){#fig:sam32-sam33-composition-part3}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM34--SAM38.](../figs/submap/case_study/SAM34-SAM35-SAM36-SAM37-SAM38-composition-part1.png){#fig:sam34-sam35-SAM36-SAM37-SAM38-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM39--SAM40.](../figs/submap/case_study/SAM39-SAM40-composition-part2.png){#fig:sam39-sam40-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM41--SAM46.](../figs/submap/case_study/SAM41-SAM42-SAM43-SAM44-SAM45-SAM46-composition.png){#fig:sam41-sam42-sam43-sam44-sam45-sam46-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM47--SAM51.](../figs/submap/case_study/SAM47-SAM48-SAM49-SAM50-SAM51-composition.png){#fig:sam47-sam48-sam49-sam50-sam51-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM52--SAM56.](../figs/submap/case_study/SAM52-SAM53-SAM54-SAM55-SAM56-composition.png){#fig:sam52-sam53-sam54-sam55-sam56-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM57--SAM61.](../figs/submap/case_study/SAM57-SAM58-SAM59-SAM60-SAM61-composition.png){#fig:sam57-sam58-sam59-sam60-sam61-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SAM62--SAM66.](../figs/submap/case_study/SAM62-SAM63-SAM64-SAM65-SAM66-composition.png){#fig:sam62-sam63-sam64-sam65-sam66-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA01--SEA05.](../figs/submap/case_study/SEA01-SEA02-SEA03-SEA04-SEA05-composition-part1.png){#fig:sea01-sea02-sea03-sea04-sea05-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA06--SEA07.](../figs/submap/case_study/SEA06-SEA07-composition-part2.png){#fig:sea06-sea07-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA08--SEA12.](../figs/submap/case_study/SEA08-SEA09-SEA10-SEA11-SEA12-composition.png){#fig:sea08-sea09-sea10-sea11-sea12-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring Java transects SEA13--SEA17.](../figs/submap/case_study/SEA13-SEA14-SEA15-SEA16-SEA17-composition.png){#fig:sea13-sea14-sea15-sea16-sea17-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA18--SEA22.](../figs/submap/case_study/SEA18-SEA19-SEA20-SEA21-SEA22-composition-part1.png){#fig:sea18-sea19-sea20-sea21-sea22-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA23--SEA27.](../figs/submap/case_study/SEA23-SEA24-SEA25-SEA26-SEA27-composition-part1.png){#fig:sea23-sea24-sea25-sea26-sea27-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA28--SEA34.](../figs/submap/case_study/SEA28-SEA31-SEA32-SEA33-SEA34-composition-part2.png){#fig:sea28-sea31-sea32-sea33-sea34-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA29--SEA30.](../figs/submap/case_study/SEA29-SEA30-composition-part2.png){#fig:sea29-sea30-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA35--SEA39.](../figs/submap/case_study/SEA35-SEA36-SEA37-SEA38-SEA39-composition-part1.png){#fig:sea35-sea36-sea37-sea38-sea39-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA40--SEA44.](../figs/submap/case_study/SEA40-SEA41-SEA42-SEA43-SEA44-composition-part2.png){#fig:sea40-sea41-sea42-sea43-sea44-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA45--SEA49.](../figs/submap/case_study/SEA45-SEA46-SEA47-SEA48-SEA49-composition-part3.png){#fig:sea45-sea46-sea47-sea48-sea49-composition-part3}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA50--SEA54.](../figs/submap/case_study/SEA50-SEA51-SEA52-SEA53-SEA54-composition-part4.png){#fig:sea50-sea51-sea52-sea53-sea54-composition-part4}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA55--SEA60.](../figs/submap/case_study/SEA55-SEA56-SEA57-SEA58-SEA59-SEA60-composition.png){#fig:sea55-sea56-sea57-sea58-sea59-sea60-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA61--SEA66.](../figs/submap/case_study/SEA61-SEA62-SEA63-SEA64-SEA65-SEA66-composition.png){#fig:sea61-sea62-sea63-sea64-sea65-sea66-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA67--SEA72.](../figs/submap/case_study/SEA67-SEA68-SEA69-SEA70-SEA71-SEA72-composition.png){#fig:sea67-sea68-sea69-sea70-sea71-sea72-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA73--SEA77.](../figs/submap/case_study/SEA73-SEA74-SEA75-SEA76-SEA77-composition.png){#fig:sea73-sea74-sea75-sea76-sea77-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA78--SEA82.](../figs/submap/case_study/SEA78-SEA79-SEA80-SEA81-SEA82-composition-part1.png){#fig:sea78-sea79-sea80-sea81-sea82-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA83--SEA86.](../figs/submap/case_study/SEA83-SEA84-SEA85-SEA86-composition-part2.png){#fig:sea83-sea84-sea85-sea86-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SEA87--SEA90.](../figs/submap/case_study/SEA87-SEA88-SEA89-SEA90-composition.png){#fig:sea87-sea88-sea89-sea90-composition}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SWP01--SWP05.](../figs/submap/case_study/SWP01-SWP02-SWP03-SWP04-SWP05-composition-part1.png){#fig:swp01-swp02-swp03-swp04-swp05-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SWP06--SWP09.](../figs/submap/case_study/SWP06-SWP07-SWP08-SWP09-composition-part2.png){#fig:swp06-swp07-swp08-swp09-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SWP10--SWP14.](../figs/submap/case_study/SWP10-SWP11-SWP12-SWP13-SWP14-composition-part1.png){#fig:swp10-swp11-swp12-swp13-swp14-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SWP15--SWP16.](../figs/submap/case_study/SWP15-SWP16-composition-part2.png){#fig:swp15-swp16-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SWP17--SWP21.](../figs/submap/case_study/SWP17-SWP18-SWP19-SWP20-SWP21-composition-part1.png){#fig:swp17-swp18-swp19-swp20-swp21-composition-part1}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SWP22--SWP26.](../figs/submap/case_study/SWP22-SWP23-SWP24-SWP25-SWP26-composition-part2.png){#fig:swp22-swp23-swp24-swp25-swp26-composition-part2}

![Cross-correlation function (CCF) analysis (top row) and surface heat flow profiles (bottom row) for neighboring transects SWP27--SWP32.](../figs/submap/case_study/SWP27-SWP28-SWP29-SWP30-SWP31-SWP32-composition.png){#fig:swp27-swp28-swp29-swp30-swp31-swp32-composition}

\clearpage

# Supplementary Data (Tables S1--S8) {#sec:tables}

|Label    |Description                                               |
|:--------|:---------------------------------------------------------|
|Table S1 |IHFC 2024 observation summary statistics per transect set |
|Table S2 |IHFC 2024 observation summary statistics per transect     |
|Table S3 |Similarity--Krige estimate differences                    |
|Table S4 |Similarity--Krige uncertainty differences                 |
|Table S5 |Similarity--Observation residuals                         |
|Table S6 |Krige--Observation residuals                              |
|Table S7 |Global adjacent-profile cross-correlations                |
|Table S8 |Per-transect method-pair cross-correlations               |
|Table S9 |Optimized Kriging parameters                              |

\clearpage

|SubMap Set | Region|    n|  Min|   Max| Median|   IQR|  Mean| $\sigma$|
|:----------|------:|----:|----:|-----:|------:|-----:|-----:|--------:|
|NPA_SET1   |    NPA| 1170|  6.9| 245.0|   80.0|  49.0|  85.1|     41.2|
|NPA_SET2   |    NPA|  624|  6.0| 245.0|   72.0|  43.6|  77.5|     40.3|
|NPA_SET3   |    NPA|  339|  6.0| 230.0|   70.0|  37.0|  78.8|     40.3|
|NPA_SET4   |    NPA|  150|  5.0| 223.0|   71.0|  39.1|  81.8|     44.3|
|NPA_SET5   |    NPA|   94|  5.0| 195.0|   71.5|  30.1|  72.4|     38.0|
|NPA_SET6   |    NPA|   43| 17.0| 129.0|   76.0|  27.5|  71.8|     28.4|
|NPA_SET7   |    NPA|   79| 17.0| 140.0|   70.0|  29.0|  74.9|     24.5|
|NPA_SET8   |    NPA| 5972|  0.6| 250.0|  100.0|  86.8| 112.2|     58.6|
|SAM_SET1   |    SAM|  663|  3.7| 250.0|   57.0|  32.0|  66.2|     45.1|
|SAM_SET2   |    SAM|  855|  3.0| 248.9|   45.9|  63.6|  60.7|     47.2|
|SAM_SET3   |    SAM|  152| 18.0| 219.0|   51.0|  32.4|  60.7|     39.5|
|SAM_SET4   |    SAM|  631|  3.0| 219.0|   53.0|  29.2|  61.9|     34.3|
|SAM_SET5   |    SAM|  612|  5.0| 250.0|  192.0| 113.7| 164.3|     72.1|
|SAM_SET6   |    SAM|   94|  6.5| 235.0|   49.0|  27.5|  61.6|     46.7|
|SAM_SET7   |    SAM|  128|  5.0| 237.0|   61.5|  60.5|  80.8|     58.1|
|SAM_SET8   |    SAM|  169|  8.0| 246.0|   60.0|  56.2|  77.4|     53.2|
|SAM_SET9   |    SAM|  204|  7.0| 246.0|   74.2|  77.5|  88.3|     56.4|
|SAM_SET10  |    SAM|   73|  2.0| 234.3|   96.7|  60.0| 106.5|     53.7|
|SEA_SET1   |    SEA|  489| 10.9| 221.0|   26.1|  28.2|  41.0|     33.0|
|SEA_SET2   |    SEA|  519| 12.5| 247.0|  103.2|  53.3|  98.9|     41.4|
|SEA_SET3   |    SEA|  230| 16.0| 199.0|   76.5|  35.2|  75.3|     28.1|
|SEA_SET4   |    SEA|  126|  1.0| 154.0|   59.2|  25.1|  58.1|     27.6|
|SEA_SET5   |    SEA|  235|  3.5| 229.3|   62.8|  37.0|  72.5|     38.1|
|SEA_SET6   |    SEA|  505|  1.0| 229.3|   59.8|  30.9|  64.8|     33.4|
|SEA_SET7   |    SEA|  439|  3.0| 231.0|   73.0|  36.7|  78.1|     39.3|
|SEA_SET8   |    SEA|  620|  3.0| 244.0|   68.3|  43.6|  80.7|     48.1|
|SEA_SET9   |    SEA| 1571| 11.0| 245.0|   82.0|  47.0|  88.1|     41.2|
|SEA_SET10  |    SEA|  778|  4.0| 242.0|   72.7|  63.0|  84.4|     50.1|
|SEA_SET11  |    SEA|  169|  1.0| 230.0|   52.0|  65.0|  67.0|     54.0|
|SEA_SET12  |    SEA|   54| 10.0| 209.0|   76.0|  44.0|  77.7|     40.7|
|SWP_SET1   |    SWP|  188|  2.0| 242.0|   51.7|  45.5|  62.1|     47.0|
|SWP_SET2   |    SWP|  101|  2.0| 234.0|   72.0|  80.5|  87.2|     60.8|
|SWP_SET3   |    SWP|  169|  5.0| 217.7|   41.0|  33.2|  51.4|     38.1|
|SWP_SET4   |    SWP|  390|  3.0| 249.0|   65.0|  52.0|  81.0|     48.8|

Table: Summary statistics of filtered and deduplicated IHFC 2024 heat flow observations (mW m$^{-2}$) distributed across the 34 SubMap transect spatial domains. NPA: North Pacific--Aleutians; SAM: South America; SEA: Southeast Asia; SWP: Southwest Pacific. {#tbl:global-obs-summary-table}

\clearpage

|SubMap Set | Transect ID|     Trench Name| Region|    n|   Min|   Max| Median|   IQR|  Mean| $\sigma$|
|:----------|-----------:|---------------:|------:|----:|-----:|-----:|------:|-----:|-----:|--------:|
|NPA_SET1   |       NPA01|           Japan|    NPA|   83|  20.1| 186.3|   81.4|  46.2|  77.2|     33.3|
|NPA_SET1   |       NPA02|           Japan|    NPA|  112|  36.0| 242.0|   93.5|  34.0|  99.1|     37.0|
|NPA_SET1   |       NPA03|           Japan|    NPA|  223|  20.0| 229.0|   76.0|  42.5|  85.6|     42.1|
|NPA_SET1   |       NPA04|           Japan|    NPA|  134|  16.5| 245.0|   86.2|  50.1|  83.9|     35.8|
|NPA_SET2   |       NPA05|          Kurils|    NPA|   78|  24.0| 226.0|   76.5|  59.2|  75.4|     39.7|
|NPA_SET2   |       NPA06|          Kurils|    NPA|   52|  36.0| 226.0|   94.0|  59.8|  96.4|     44.3|
|NPA_SET2   |       NPA07|          Kurils|    NPA|   53|  22.0| 160.0|   83.0|  44.0|  84.1|     32.9|
|NPA_SET2   |       NPA08|          Kurils|    NPA|   31|  25.0| 188.0|   80.0|  46.4|  81.9|     41.0|
|NPA_SET2   |       NPA09|          Kurils|    NPA|   31|  16.8| 188.0|   67.0|  38.5|  74.1|     38.3|
|NPA_SET2   |       NPA10|          Kurils|    NPA|   39|   6.0| 110.0|   59.0|  33.2|  61.3|     24.6|
|NPA_SET3   |       NPA11|          Kurils|    NPA|   41|  22.2| 230.0|   80.0|  38.1|  86.9|     42.4|
|NPA_SET3   |       NPA12|       Kamchatka|    NPA|   25|  29.0| 173.0|   69.0|  20.3|  73.8|     31.5|
|NPA_SET3   |       NPA13|       Kamchatka|    NPA|   26|  30.0| 154.0|   61.5|  19.9|  64.3|     22.4|
|NPA_SET3   |       NPA14|       Kamchatka|    NPA|   19|  33.5| 147.0|   72.0|  54.6|  75.4|     35.1|
|NPA_SET4   |       NPA15|       Aleutians|    NPA|   28|  16.0| 223.0|  108.0|  45.7| 112.8|     46.6|
|NPA_SET4   |       NPA16|       Aleutians|    NPA|   21|  28.0| 171.0|   58.0|  21.6|  70.5|     36.2|
|NPA_SET4   |       NPA17|       Aleutians|    NPA|   14|  46.0|  87.0|   69.5|  20.8|  67.1|     15.3|
|NPA_SET4   |       NPA18|       Aleutians|    NPA|   11|   5.0|  77.3|   62.0|  19.0|  58.7|     20.8|
|NPA_SET4   |       NPA19|       Aleutians|    NPA|    8|  45.0|  73.1|   67.0|  19.1|  62.3|     12.2|
|NPA_SET4   |       NPA20|       Aleutians|    NPA|   12|  45.0|  87.0|   73.0|   9.8|  71.0|     13.8|
|NPA_SET4   |       NPA21|       Aleutians|    NPA|    6|  56.0| 161.0|   77.1|  25.5|  92.5|     37.6|
|NPA_SET4   |       NPA22|       Aleutians|    NPA|   12|  34.0| 109.0|   61.8|  32.0|  64.9|     24.7|
|NPA_SET5   |       NPA23|       Aleutians|    NPA|    8|  43.9|  77.2|   62.8|  11.3|  62.4|     11.1|
|NPA_SET5   |       NPA24|       Aleutians|    NPA|   12|  35.9| 195.0|   66.0|  75.2|  95.9|     62.4|
|NPA_SET5   |       NPA25|       Aleutians|    NPA|    9|  38.0|  81.0|   49.0|  22.1|  54.7|     16.6|
|NPA_SET5   |       NPA26|       Aleutians|    NPA|    8|   7.5| 132.0|   66.0|  52.9|  67.6|     48.2|
|NPA_SET5   |       NPA27|       Aleutians|    NPA|   10|  53.0| 195.0|   86.5|  20.7|  95.8|     37.7|
|NPA_SET5   |       NPA28|       Aleutians|    NPA|    6|   5.8| 195.0|   76.0|   4.2|  84.2|     61.2|
|NPA_SET5   |       NPA29|       Aleutians|    NPA|    7|  30.5| 113.0|   44.0|  59.0|  64.8|     36.7|
|NPA_SET5   |       NPA30|       Aleutians|    NPA|    4|  78.0| 122.0|  100.2|  43.8| 100.1|     25.3|
|NPA_SET6   |       NPA31|          Alaska|    NPA|    1|  17.0|  17.0|   17.0|   0.0|  17.0|         |
|NPA_SET6   |       NPA32|          Alaska|    NPA|    7|  26.0|  81.0|   67.0|  17.0|  62.6|     18.4|
|NPA_SET6   |       NPA33|          Alaska|    NPA|    8|  17.0| 129.0|   75.0|  22.8|  70.6|     31.8|
|NPA_SET6   |       NPA34|          Alaska|    NPA|    5|  17.0| 113.0|   77.0|  14.8|  70.8|     34.8|
|NPA_SET6   |       NPA35|          Alaska|    NPA|    5|  66.0|  94.0|   87.0|  16.0|  83.4|     11.8|
|NPA_SET6   |       NPA36|          Alaska|    NPA|    3|  24.7|  74.0|   71.0|  24.6|  56.6|     27.6|
|NPA_SET7   |       NPA37|          Alaska|    NPA|    3|  70.0| 106.0|   87.0|  18.0|  87.7|     18.0|
|NPA_SET7   |       NPA38|          Alaska|    NPA|    5|  46.0| 106.0|   70.0|  30.0|  73.2|     23.9|
|NPA_SET7   |       NPA39|          Alaska|    NPA|    5|  50.0| 106.0|   61.0|   2.1|  67.4|     22.0|
|NPA_SET7   |       NPA40|          Alaska|    NPA|   12|  59.0| 140.0|   70.0|  27.7|  83.0|     28.9|
|NPA_SET7   |       NPA41|          Alaska|    NPA|   11|  60.0| 140.0|   71.0|  30.0|  87.7|     29.2|
|NPA_SET7   |       NPA42|          Alaska|    NPA|   21|  57.0| 140.0|   73.0|  25.0|  82.4|     24.4|
|NPA_SET8   |       NPA43|        Cascadia|    NPA|  127|  26.0| 236.3|  106.0|  73.0| 111.5|     53.0|
|NPA_SET8   |       NPA44|        Cascadia|    NPA| 1077|   3.0| 250.0|  157.0| 102.0| 151.9|     60.5|
|NPA_SET8   |       NPA45|        Cascadia|    NPA|  293|   9.0| 244.0|   60.0|  42.0|  72.6|     42.1|
|NPA_SET8   |       NPA46|        Cascadia|    NPA|  790|   6.0| 249.0|   91.0|  54.0|  97.1|     44.6|
|NPA_SET8   |       NPA47|        Cascadia|    NPA|  503|   1.0| 246.9|   83.0|  76.5|  93.0|     55.0|
|SAM_SET1   |       SAM01| Central America|    SAM|   23|  35.0| 249.3|   59.0|  73.0|  95.0|     62.8|
|SAM_SET1   |       SAM02| Central America|    SAM|   44|  13.0| 223.0|   57.0|  35.2|  72.4|     49.3|
|SAM_SET1   |       SAM03| Central America|    SAM|   60|  34.0| 248.1|   61.0| 117.2| 101.7|     66.6|
|SAM_SET1   |       SAM04| Central America|    SAM|   68|  13.0| 245.1|   60.5|  19.5|  73.1|     51.0|
|SAM_SET1   |       SAM05| Central America|    SAM|   14|  27.0| 106.0|   53.0|  50.2|  59.9|     29.5|
|SAM_SET1   |       SAM06| Central America|    SAM|   24|  14.0| 158.0|   36.5|  27.7|  51.1|     46.0|
|SAM_SET1   |       SAM07| Central America|    SAM|   43|  15.0| 147.0|   45.0|  16.0|  46.9|     22.5|
|SAM_SET1   |       SAM08| Central America|    SAM|    6|  13.0|  54.0|   48.0|  25.4|  39.2|     18.0|
|SAM_SET2   |       SAM09| Central America|    SAM|   64|   3.7|  73.6|   25.4|  36.6|  31.4|     20.4|
|SAM_SET2   |       SAM10| Central America|    SAM|   45|   4.0|  95.0|   26.0|  12.0|  27.2|     15.5|
|SAM_SET2   |       SAM11| Central America|    SAM|  387|   3.0| 244.0|   37.0|  79.7|  59.4|     46.5|
|SAM_SET2   |       SAM12| Central America|    SAM|   96|   7.9| 248.9|   60.7|  37.4|  78.8|     47.9|
|SAM_SET2   |       SAM13|          Panama|    SAM|    7|  97.6| 176.0|  156.0|  15.0| 151.9|     26.4|
|SAM_SET2   |       SAM14|          Panama|    SAM|   20|  41.0| 176.0|   76.0|  58.6|  94.0|     43.9|
|SAM_SET2   |       SAM15|          Panama|    SAM|   24|  17.0| 165.0|   73.4|  75.7|  82.5|     45.7|
|SAM_SET3   |       SAM16|       Venezuela|    SAM|    5|  53.6|  90.0|   63.0|   3.0|  66.7|     13.7|
|SAM_SET3   |       SAM17|       Venezuela|    SAM|   11|  18.0|  90.0|   51.0|  20.7|  51.4|     21.5|
|SAM_SET3   |       SAM18|       Venezuela|    SAM|    4|  52.7| 129.0|   90.2|  72.8|  90.5|     42.7|
|SAM_SET3   |       SAM19|       Venezuela|    SAM|   18|  19.0| 129.0|   46.1|  19.8|  49.4|     25.9|
|SAM_SET3   |       SAM20|       Venezuela|    SAM|    6|  44.0| 107.0|   49.5|   7.8|  58.0|     24.3|
|SAM_SET3   |       SAM21|       Venezuela|    SAM|   13|  22.0| 200.0|   81.0|  78.9|  99.9|     59.0|
|SAM_SET4   |       SAM22|         Muertos|    SAM|    9|  15.0|  53.0|   44.0|  31.0|  35.0|     17.2|
|SAM_SET4   |       SAM23|         Muertos|    SAM|   53|  37.3|  74.0|   47.7|   3.4|  48.6|      6.2|
|SAM_SET4   |       SAM24|         Muertos|    SAM|   13|  29.0|  72.0|   54.0|  17.0|  53.7|     12.5|
|SAM_SET4   |       SAM25|      Hispaniola|    SAM|   13|  15.0|  67.0|   50.0|  36.7|  40.7|     18.9|
|SAM_SET4   |       SAM26|     Puerto-Rico|    SAM|    4|  25.0|  74.0|   54.5|  14.5|  52.0|     20.2|
|SAM_SET4   |       SAM27|     Puerto-Rico|    SAM|   12|  29.0|  72.0|   52.6|  16.8|  53.8|     12.8|
|SAM_SET4   |       SAM28|        Antilles|    SAM|   16|  43.0| 219.0|   87.0|  34.7|  98.8|     51.9|
|SAM_SET4   |       SAM29|        Antilles|    SAM|   25|  22.5| 219.0|   57.0|  35.5|  64.7|     50.1|
|SAM_SET4   |       SAM30|        Antilles|    SAM|   53|  11.0| 219.0|   95.0|  38.1|  96.8|     41.8|
|SAM_SET4   |       SAM31|        Antilles|    SAM|  129|  29.0| 219.0|   70.0|  24.6|  73.6|     28.7|
|SAM_SET4   |       SAM32|        Antilles|    SAM|  175|   3.0| 214.0|   45.0|  20.0|  48.8|     32.0|
|SAM_SET4   |       SAM33|        Antilles|    SAM|   27|  23.0| 118.0|   57.0|  16.0|  58.2|     19.2|
|SAM_SET5   |       SAM34|          Andean|    SAM|    9|  23.2| 155.0|   81.5|  78.9|  87.3|     50.5|
|SAM_SET5   |       SAM35|          Andean|    SAM|    4|  79.0| 147.0|   98.2|  17.3| 105.6|     29.0|
|SAM_SET5   |       SAM36|          Andean|    SAM|   21|  65.0| 225.0|   72.0| 106.0| 107.6|     58.8|
|SAM_SET5   |       SAM37|          Andean|    SAM|   11|  29.0| 219.0|   70.0|  86.0|  98.9|     61.3|
|SAM_SET5   |       SAM38|          Andean|    SAM|    3|  50.0|  53.0|   52.0|   1.5|  51.7|      1.5|
|SAM_SET5   |       SAM39|          Andean|    SAM|    5|  28.0|  82.0|   37.0|  15.0|  46.4|     21.4|
|SAM_SET5   |       SAM40|          Andean|    SAM|    1|  55.0|  55.0|   55.0|   0.0|  55.0|         |
|SAM_SET6   |       SAM41|          Andean|    SAM|    4|  39.0|  55.0|   49.5|  12.2|  48.2|      8.1|
|SAM_SET6   |       SAM42|          Andean|    SAM|   12|  29.0|  69.0|   37.5|  11.0|  40.7|     10.7|
|SAM_SET6   |       SAM43|          Andean|    SAM|    3|  49.0|  77.0|   49.0|  14.0|  58.3|     16.2|
|SAM_SET6   |       SAM44|          Andean|    SAM|    8|   6.5| 196.0|   49.9|  64.3|  62.0|     65.2|
|SAM_SET6   |       SAM45|          Andean|    SAM|    2|  12.0|  27.0|   19.5|   7.5|  19.5|     10.6|
|SAM_SET6   |       SAM46|          Andean|    SAM|    2|  44.0|  64.0|   54.0|  10.0|  54.0|     14.1|
|SAM_SET7   |       SAM47|          Andean|    SAM|    4|  65.0| 228.0|  176.0|  51.2| 161.2|     68.9|
|SAM_SET7   |       SAM48|          Andean|    SAM|   10|   5.0| 156.0|   46.5|  40.2|  64.5|     51.0|
|SAM_SET7   |       SAM49|          Andean|    SAM|    8|  10.0|  47.0|   27.5|  12.3|  28.0|     11.6|
|SAM_SET7   |       SAM50|          Andean|    SAM|    4| 135.0| 215.0|  172.5|  68.8| 173.8|     42.1|
|SAM_SET7   |       SAM51|          Andean|    SAM|    9|  10.0| 105.0|   77.0|  45.0|  59.2|     36.8|
|SAM_SET8   |       SAM52|          Andean|    SAM|   10|  21.0| 194.0|  162.0| 105.0| 127.3|     68.2|
|SAM_SET8   |       SAM53|          Andean|    SAM|    2|  72.0| 166.0|  119.0|  47.0| 119.0|     66.5|
|SAM_SET8   |       SAM54|          Andean|    SAM|    5|  20.9|  70.0|   31.0|  37.6|  41.9|     22.7|
|SAM_SET8   |       SAM55|          Andean|    SAM|    9|  57.0| 212.0|  102.0|  93.0| 123.3|     57.0|
|SAM_SET8   |       SAM56|          Andean|    SAM|   52|  30.0| 190.0|   53.0|  33.0|  64.5|     34.0|
|SAM_SET9   |       SAM57|          Andean|    SAM|   18|  46.0| 215.0|   87.5|  90.8| 111.5|     56.4|
|SAM_SET9   |       SAM58|          Andean|    SAM|    9|  23.0| 168.0|   60.0|  86.0|  82.9|     51.3|
|SAM_SET9   |       SAM59|          Andean|    SAM|    9|  27.0| 218.0|   39.0|  14.0|  74.8|     76.4|
|SAM_SET9   |       SAM60|          Andean|    SAM|   32|   7.0| 200.6|   75.3|  62.4|  85.3|     48.4|
|SAM_SET9   |       SAM61|          Andean|    SAM|   27|  27.0| 234.3|   77.4|  64.0| 105.0|     64.7|
|SAM_SET10  |       SAM62|          Andean|    SAM|   14|  79.0| 148.0|  122.5|  34.5| 121.7|     24.4|
|SAM_SET10  |       SAM63|          Andean|    SAM|    0|      |      |       |      |      |         |
|SAM_SET10  |       SAM64|          Andean|    SAM|    0|      |      |       |      |      |         |
|SAM_SET10  |       SAM65|          Andean|    SAM|    1|  96.7|  96.7|   96.7|   0.0|  96.7|         |
|SAM_SET10  |       SAM66|          Andean|    SAM|    1|  33.9|  33.9|   33.9|   0.0|  33.9|         |
|SEA_SET1   |       SEA01|         Andaman|    SEA|    1|  81.0|  81.0|   81.0|   0.0|  81.0|         |
|SEA_SET1   |       SEA02|         Andaman|    SEA|    3|  38.0| 100.0|   84.4|  31.0|  74.1|     32.2|
|SEA_SET1   |       SEA03|         Andaman|    SEA|   41|  10.9| 221.0|   22.4|  19.1|  27.4|     32.5|
|SEA_SET1   |       SEA04|         Andaman|    SEA|    3|  43.0|  84.4|   43.4|  20.7|  56.9|     23.8|
|SEA_SET1   |       SEA05|         Andaman|    SEA|    0|      |      |       |      |      |         |
|SEA_SET1   |       SEA06|         Sumatra|    SEA|    5|  61.0| 105.0|   73.0|  29.5|  80.3|     19.2|
|SEA_SET1   |       SEA07|         Sumatra|    SEA|   13|  29.0| 106.6|   78.7|  25.5|  68.8|     21.8|
|SEA_SET2   |       SEA08|         Sumatra|    SEA|   67|  28.0| 226.0|  108.0|  68.6|  97.3|     40.4|
|SEA_SET2   |       SEA09|         Sumatra|    SEA|  127|  37.6| 237.8|  128.0|  34.6| 132.7|     34.5|
|SEA_SET2   |       SEA10|         Sumatra|    SEA|   52|  12.5| 167.2|  108.0|  30.9|  95.9|     34.1|
|SEA_SET2   |       SEA11|         Sumatra|    SEA|   36|  28.0| 128.3|   96.8|  52.3|  81.7|     31.3|
|SEA_SET2   |       SEA12|         Sumatra|    SEA|   34|  24.0| 199.0|  101.2|  24.2|  92.7|     33.9|
|SEA_SET3   |       SEA13|            Java|    SEA|   31|  58.1| 168.0|   77.3|  20.3|  85.5|     25.6|
|SEA_SET3   |       SEA14|            Java|    SEA|   42|  45.6| 149.6|   79.6|  18.3|  80.8|     20.9|
|SEA_SET3   |       SEA15|            Java|    SEA|    8|  20.0|  93.2|   48.9|  10.8|  53.0|     21.1|
|SEA_SET3   |       SEA16|            Java|    SEA|    7|  20.9| 115.0|   91.5|  18.3|  83.5|     30.5|
|SEA_SET3   |       SEA17|            Java|    SEA|   25|  51.8| 117.9|   78.2|  21.3|  78.1|     16.9|
|SEA_SET4   |       SEA18|            Java|    SEA|    9|  16.0|  92.8|   74.8|  43.0|  61.4|     27.7|
|SEA_SET4   |       SEA19|            Java|    SEA|    3|  39.0|  76.0|   46.0|  18.5|  53.7|     19.7|
|SEA_SET4   |       SEA20|           Timor|    SEA|   10|   3.5|  66.0|   17.3|  39.9|  25.6|     22.8|
|SEA_SET4   |       SEA21|           Timor|    SEA|    2|  64.0|  71.0|   67.5|   3.5|  67.5|      4.9|
|SEA_SET4   |       SEA22|           Timor|    SEA|    6|  50.0|  63.0|   53.5|   5.5|  55.2|      4.9|
|SEA_SET4   |       SEA29|          Flores|    SEA|    3|  39.0|  70.0|   48.0|  15.5|  52.3|     15.9|
|SEA_SET4   |       SEA30|          Flores|    SEA|    4|   6.9|  67.0|   56.0|  21.0|  46.5|     27.1|
|SEA_SET5   |       SEA23|           Timor|    SEA|    9|  48.0| 103.0|   58.0|  13.0|  65.7|     21.8|
|SEA_SET5   |       SEA24|           Timor|    SEA|    8|  53.0| 154.0|   71.5|  64.0|  94.2|     42.3|
|SEA_SET5   |       SEA25|           Timor|    SEA|    6|   4.0| 148.0|  113.0|  47.8|  96.8|     52.3|
|SEA_SET5   |       SEA26|        Tanimbar|    SEA|   16|  46.0| 227.0|   69.4|  56.3| 101.4|     60.2|
|SEA_SET5   |       SEA27|        Tanimbar|    SEA|   20|   4.0| 175.0|   67.0|  64.6|  80.8|     46.1|
|SEA_SET5   |       SEA28|             Aru|    SEA|   18|   4.0| 175.0|   70.3|  71.8|  79.9|     47.2|
|SEA_SET5   |       SEA31|           Wetar|    SEA|    9|  33.0|  94.1|   53.0|  12.0|  56.9|     17.4|
|SEA_SET5   |       SEA32|           Seram|    SEA|   46|  20.9| 169.5|   58.6|  21.5|  66.4|     35.0|
|SEA_SET5   |       SEA33|           Seram|    SEA|   32|  16.7| 122.1|   57.0|  23.0|  55.8|     23.6|
|SEA_SET5   |       SEA34|           Seram|    SEA|    7|  46.0| 107.0|   71.0|  12.7|  73.8|     18.7|
|SEA_SET6   |       SEA35|       Halmahera|    SEA|   23|  48.9| 229.3|   80.8|  24.8|  87.9|     34.5|
|SEA_SET6   |       SEA36|       Halmahera|    SEA|    6|  52.7|  95.4|   59.2|  25.4|  69.3|     18.8|
|SEA_SET6   |       SEA37|  North Sulawesi|    SEA|   14|  12.5| 154.0|   67.1|  31.0|  72.5|     42.6|
|SEA_SET6   |       SEA38|  North Sulawesi|    SEA|   24|   1.0|  94.7|   11.7|  43.1|  25.6|     26.2|
|SEA_SET6   |       SEA39|         Sangihe|    SEA|   10|  11.2| 223.4|   68.2|  37.6|  88.1|     57.8|
|SEA_SET6   |       SEA40|         Sangihe|    SEA|    9|  41.0| 101.0|   71.0|   8.8|  69.9|     21.3|
|SEA_SET6   |       SEA41|         Sangihe|    SEA|    9|  32.0| 101.0|   67.0|  23.6|  68.9|     25.4|
|SEA_SET6   |       SEA42|     Philippines|    SEA|    3|  31.0|  46.0|   37.0|   7.5|  38.0|      7.5|
|SEA_SET6   |       SEA43|     Philippines|    SEA|    7|  17.2|  71.2|   69.3|   8.7|  60.7|     19.6|
|SEA_SET6   |       SEA44|     Philippines|    SEA|    8|  29.3| 109.0|   80.8|  47.5|  71.5|     28.6|
|SEA_SET6   |       SEA45|     Philippines|    SEA|   11|  24.7| 109.0|   79.6|  39.5|  71.4|     27.6|
|SEA_SET6   |       SEA46|     Philippines|    SEA|    7|  17.0|  62.0|   42.7|  13.3|  40.2|     14.2|
|SEA_SET6   |       SEA47|     Philippines|    SEA|    4|  45.0|  83.0|   62.5|  10.2|  63.2|     15.5|
|SEA_SET6   |       SEA48|     Philippines|    SEA|    5|  17.0|  82.1|   45.0|  15.9|  48.0|     23.8|
|SEA_SET6   |       SEA49|        Cotobato|    SEA|   10|  10.0|  71.2|   62.7|  37.8|  49.0|     25.3|
|SEA_SET6   |       SEA50|        Cotobato|    SEA|    9|  29.3|  82.1|   67.0|   3.9|  65.1|     14.6|
|SEA_SET6   |       SEA51|            Sulu|    SEA|    2|  67.8|  69.1|   68.5|   0.6|  68.5|      0.9|
|SEA_SET6   |       SEA52|            Sulu|    SEA|    9|  17.0|  96.0|   67.0|  19.0|  70.8|     23.6|
|SEA_SET6   |       SEA53|            Sulu|    SEA|   11|  11.2| 223.4|   69.3|  33.1|  77.8|     57.0|
|SEA_SET6   |       SEA54|          Manila|    SEA|    2|  40.0|  55.9|   48.0|   7.9|  48.0|     11.2|
|SEA_SET7   |       SEA55|          Manila|    SEA|    7|  13.0| 117.0|   55.9|  73.0|  66.8|     42.8|
|SEA_SET7   |       SEA56|          Manila|    SEA|   11|  22.9| 152.0|   89.3|  48.0|  74.8|     37.7|
|SEA_SET7   |       SEA57|          Manila|    SEA|   19|  49.0| 146.0|   81.0|  16.0|  82.7|     20.1|
|SEA_SET7   |       SEA58|          Manila|    SEA|   26|  16.0| 147.0|   77.0|  25.9|  75.4|     31.1|
|SEA_SET7   |       SEA59|          Manila|    SEA|   48|   3.0| 219.3|   61.9|  29.2|  71.7|     44.3|
|SEA_SET7   |       SEA60|          Manila|    SEA|   51|  16.0| 209.0|   74.0|  39.9|  76.9|     39.3|
|SEA_SET8   |       SEA61|         Ryukyus|    SEA|   41|   8.0| 182.1|   40.0|  27.5|  46.7|     34.4|
|SEA_SET8   |       SEA62|         Ryukyus|    SEA|   21|   8.8| 182.1|   68.0|  31.1|  66.0|     37.7|
|SEA_SET8   |       SEA63|         Ryukyus|    SEA|   22|  28.0| 227.0|   86.0|  77.5| 102.8|     49.9|
|SEA_SET8   |       SEA64|         Ryukyus|    SEA|  118|   9.0| 244.0|  101.0|  91.8| 110.2|     65.4|
|SEA_SET8   |       SEA65|         Ryukyus|    SEA|   20|  16.0| 168.0|   70.0|  54.6|  83.7|     42.2|
|SEA_SET8   |       SEA66|         Ryukyus|    SEA|   18|  29.0| 120.0|   55.5|  28.6|  57.5|     22.0|
|SEA_SET9   |       SEA67|         Ryukyus|    SEA|  125|  32.0| 228.0|   65.0|  19.7|  68.5|     24.4|
|SEA_SET9   |       SEA68|          Nankai|    SEA|   83|  25.0| 166.0|   90.0|  44.0|  88.7|     30.3|
|SEA_SET9   |       SEA69|          Nankai|    SEA|  198|  33.0| 242.0|  106.5|  67.8| 116.4|     49.0|
|SEA_SET9   |       SEA70|          Nankai|    SEA|  176|  30.0| 224.0|   90.0|  34.9|  92.3|     32.7|
|SEA_SET9   |       SEA71|          Suruga|    SEA|  197|  15.0| 232.0|   82.0|  44.0|  85.3|     37.2|
|SEA_SET9   |       SEA72|          Sagami|    SEA|  279|  15.0| 245.0|   81.0|  71.7|  95.2|     53.4|
|SEA_SET10  |       SEA73|       Izu-Bonin|    SEA|  190|  25.0| 242.0|  120.0|  61.9| 124.2|     46.4|
|SEA_SET10  |       SEA74|       Izu-Bonin|    SEA|   61|   7.0| 205.0|   50.0|  67.0|  66.9|     48.9|
|SEA_SET10  |       SEA75|       Izu-Bonin|    SEA|   57|   8.0| 206.0|   45.3|  62.0|  59.9|     42.0|
|SEA_SET10  |       SEA76|       Izu-Bonin|    SEA|   50|   7.0| 217.0|   80.8|  52.0|  78.5|     40.2|
|SEA_SET10  |       SEA77|       Izu-Bonin|    SEA|   10|   6.0| 237.0|   29.0|  68.8|  60.5|     69.9|
|SEA_SET11  |       SEA78|         Mariana|    SEA|   36|   1.0| 208.0|   48.5|  70.2|  68.0|     63.3|
|SEA_SET11  |       SEA79|         Mariana|    SEA|   13|  16.0| 130.2|   61.0|  55.2|  57.4|     33.7|
|SEA_SET11  |       SEA80|         Mariana|    SEA|   33|   2.5| 180.0|   64.0|  66.0|  66.8|     45.4|
|SEA_SET11  |       SEA81|         Mariana|    SEA|   26|  13.0| 188.0|   63.5|  38.1|  76.9|     47.5|
|SEA_SET11  |       SEA82|         Mariana|    SEA|   13|  18.0| 130.2|   64.0|  56.1|  65.3|     38.4|
|SEA_SET11  |       SEA83|         Mariana|    SEA|   10|  18.0|  91.0|   68.5|  51.8|  58.7|     29.4|
|SEA_SET11  |       SEA84|         Mariana|    SEA|    7|  19.3|  91.0|   87.1|  21.2|  73.4|     26.0|
|SEA_SET11  |       SEA85|         Mariana|    SEA|   13|   8.0|  91.0|   77.0|  55.1|  59.8|     31.5|
|SEA_SET11  |       SEA86|         Mariana|    SEA|    8|  18.0| 130.2|   70.5|  33.7|  68.1|     36.6|
|SEA_SET12  |       SEA87|             Yap|    SEA|   16|  29.0| 175.0|   81.0|  40.0|  87.6|     35.5|
|SEA_SET12  |       SEA88|             Yap|    SEA|    4|  46.1| 209.0|  117.2|  74.1| 122.4|     69.2|
|SEA_SET12  |       SEA89|             Yap|    SEA|    4|  10.6|  86.0|   44.2|  19.1|  46.2|     30.9|
|SEA_SET12  |       SEA90|           Palau|    SEA|    5|  30.0| 116.0|   44.4|  38.1|  63.3|     35.2|
|SWP_SET1   |       SWP01|     New Britain|    SWP|    6|  37.7| 193.0|   60.9|  89.1|  89.4|     65.4|
|SWP_SET1   |       SWP02|     New Britain|    SWP|    9|  12.0| 193.0|   38.0|  18.3|  51.6|     55.6|
|SWP_SET1   |       SWP03|     New Britain|    SWP|    9|   3.0| 234.0|  138.0| 122.0| 141.8|     79.0|
|SWP_SET1   |       SWP04|        Solomons|    SWP|    9|  19.0|  99.0|   38.0|  35.0|  49.9|     25.7|
|SWP_SET1   |       SWP05|        Solomons|    SWP|    9|  20.0|  58.0|   33.0|  12.3|  35.2|     14.0|
|SWP_SET1   |       SWP06|        Solomons|    SWP|   33|   9.0| 170.0|   33.7|  31.8|  44.8|     33.0|
|SWP_SET1   |       SWP07|        Solomons|    SWP|    3|  26.0|  60.6|   55.0|  17.3|  47.2|     18.6|
|SWP_SET1   |       SWP08|        Solomons|    SWP|    0|      |      |       |      |      |         |
|SWP_SET1   |       SWP09|        Solomons|    SWP|    5|  35.0|  79.0|   58.0|  19.0|  61.0|     17.7|
|SWP_SET2   |       SWP10|    New Hebrides|    SWP|   10|  30.0| 134.8|   87.6|  59.0|  86.8|     35.9|
|SWP_SET2   |       SWP11|    New Hebrides|    SWP|    9|   4.0| 234.0|   46.0|  72.4|  61.6|     72.0|
|SWP_SET2   |       SWP12|    New Hebrides|    SWP|   13|   4.0| 234.0|   49.0| 149.0|  85.1|     86.5|
|SWP_SET2   |       SWP13|    New Hebrides|    SWP|   12|  13.0| 188.0|   95.7|  89.8| 104.1|     60.7|
|SWP_SET2   |       SWP14|    New Hebrides|    SWP|    8|  13.0| 193.0|   82.5|  46.9|  97.3|     53.4|
|SWP_SET2   |       SWP15|    New Hebrides|    SWP|    7|  54.0| 141.0|   72.0|  15.7|  76.6|     29.6|
|SWP_SET2   |       SWP16|    New Hebrides|    SWP|   11|   2.0| 141.0|   60.0|  80.2|  73.4|     49.0|
|SWP_SET3   |       SWP17|           Tonga|    SWP|    9|  15.0|  72.0|   38.0|  27.0|  42.1|     18.9|
|SWP_SET3   |       SWP18|           Tonga|    SWP|    5|  18.0| 136.0|   64.0|  50.0|  77.6|     46.1|
|SWP_SET3   |       SWP19|           Tonga|    SWP|    8|  22.0| 140.0|   47.5|  25.5|  55.1|     37.0|
|SWP_SET3   |       SWP20|           Tonga|    SWP|   11|  12.0|  71.2|   47.0|  23.0|  45.2|     20.0|
|SWP_SET3   |       SWP21|           Tonga|    SWP|    6|  26.0| 100.0|   44.5|  26.2|  50.5|     27.4|
|SWP_SET3   |       SWP22|           Tonga|    SWP|    9|  13.0| 162.0|   30.0|  15.6|  45.3|     47.6|
|SWP_SET3   |       SWP23|        Kermadec|    SWP|    4|   5.0|  98.0|   51.5|  27.0|  51.5|     38.0|
|SWP_SET3   |       SWP24|        Kermadec|    SWP|    2| 137.0| 137.0|  137.0|   0.0| 137.0|      0.0|
|SWP_SET3   |       SWP25|        Kermadec|    SWP|    0|      |      |       |      |      |         |
|SWP_SET3   |       SWP26|        Kermadec|    SWP|    3|  15.0|  68.2|   68.0|  26.6|  50.4|     30.7|
|SWP_SET4   |       SWP27|        Kermadec|    SWP|    1|  24.0|  24.0|   24.0|   0.0|  24.0|         |
|SWP_SET4   |       SWP28|        Kermadec|    SWP|    0|      |      |       |      |      |         |
|SWP_SET4   |       SWP29|       Hikurangi|    SWP|  145|   8.0| 244.0|   95.0|  63.0| 100.3|     51.2|
|SWP_SET4   |       SWP30|       Hikurangi|    SWP|  113|   3.0| 249.0|   62.0|  44.1|  77.6|     53.0|
|SWP_SET4   |       SWP31|       Hikurangi|    SWP|  103|   3.0| 249.0|   63.0|  53.1|  80.0|     55.5|
|SWP_SET4   |       SWP32|       Hikurangi|    SWP|   32|  36.0|  81.0|   48.5|  16.3|  51.1|     11.3|

Table: Summary statistics of surface heat flow observations (mW m$^{-2}$) distributed across all SubMap transects. NPA: North Pacific--Aleutians; SAM: South America; SEA: Southeast Asia; SWP: Southwest Pacific. {#tbl:submap-obs-summary-table}

\clearpage

|SubMap Set | Region|    n|    Min|   Max| Median|  IQR| Mean| $\sigma$|
|:----------|------:|----:|------:|-----:|------:|----:|----:|--------:|
|NPA_SET1   |    NPA|  715| -144.8| 147.9|    2.4| 20.7|  1.0|     28.7|
|NPA_SET2   |    NPA|  913|  -95.8| 179.6|    3.2| 22.1|  3.8|     24.0|
|NPA_SET3   |    NPA|  758| -143.8| 149.4|    3.6| 31.5|  0.5|     34.3|
|NPA_SET4   |    NPA|  755| -149.3|  79.3|    6.8| 30.8| -0.5|     33.4|
|NPA_SET5   |    NPA|  912| -111.6| 128.7|    7.8| 26.7|  6.3|     25.8|
|NPA_SET6   |    NPA|  401|  -83.5|  75.2|    5.3| 37.7|  4.9|     29.6|
|NPA_SET7   |    NPA|  650|  -49.2| 112.3|    4.3| 22.1|  5.3|     21.4|
|NPA_SET8   |    NPA|  997| -164.5| 164.6|    2.1| 28.1|  9.0|     34.1|
|SAM_SET1   |    SAM|  167| -186.6| 227.6|   17.7| 56.2| 14.8|     66.9|
|SAM_SET2   |    SAM|  981|  -73.5| 214.3|   13.6| 41.8| 20.9|     35.5|
|SAM_SET3   |    SAM|  661| -110.1| 141.6|    3.0| 25.3|  0.8|     30.5|
|SAM_SET4   |    SAM|  666| -185.2| 164.9|    1.6| 21.0| -6.9|     43.2|
|SAM_SET5   |    SAM|  796| -122.0| 190.3|   18.2| 35.2| 16.2|     33.3|
|SAM_SET6   |    SAM|  348| -173.8|  99.3|   19.5| 37.9| 19.3|     32.8|
|SAM_SET7   |    SAM|  583| -149.8| 179.1|   10.9| 43.7| 12.6|     42.4|
|SAM_SET8   |    SAM|  139| -178.3| 107.7|   -0.6| 43.1| -3.6|     45.2|
|SAM_SET9   |    SAM|  327| -139.1| 176.6|    8.2| 73.1| 10.0|     57.0|
|SAM_SET10  |    SAM|  329| -100.3| 172.3|   21.4| 73.2| 29.1|     50.4|
|SEA_SET1   |    SEA|  734| -116.4|  74.9|   -0.2| 32.3| -4.0|     24.6|
|SEA_SET2   |    SEA|  528| -128.3|  78.5|    0.6| 28.9| -3.6|     28.6|
|SEA_SET3   |    SEA|  442| -143.1|  70.3|    2.1| 20.5| -0.3|     24.8|
|SEA_SET4   |    SEA|  648|  -66.4|  94.4|    5.4| 19.7|  7.0|     21.4|
|SEA_SET5   |    SEA|  502| -174.9| 117.3|    4.2| 28.5|  5.7|     31.5|
|SEA_SET6   |    SEA| 1601| -107.4| 154.2|   20.4| 36.9| 21.1|     28.9|
|SEA_SET7   |    SEA|  600| -120.9| 148.8|    6.9| 31.9|  7.6|     36.3|
|SEA_SET8   |    SEA|  859| -120.9| 173.6|    6.8| 20.8|  6.1|     30.0|
|SEA_SET9   |    SEA|  874| -110.6| 168.8|    4.4| 20.7|  3.9|     29.4|
|SEA_SET10  |    SEA|  767|  -89.0| 126.8|    2.1| 37.6|  2.0|     31.6|
|SEA_SET11  |    SEA|  359| -116.7| 121.2|   20.7| 51.0| 20.4|     45.5|
|SEA_SET12  |    SEA|  242| -151.7| 168.4|    8.1| 42.5| 11.7|     43.5|
|SWP_SET1   |    SWP|  593| -137.2| 116.8|    5.0| 30.6|  5.1|     33.2|
|SWP_SET2   |    SWP|  551| -168.3| 167.2|   12.4| 66.7| 16.6|     56.3|
|SWP_SET3   |    SWP|  518|  -89.1| 146.0|   13.7| 44.4| 20.4|     36.6|
|SWP_SET4   |    SWP|  836| -101.0| 141.1|   13.4| 33.3| 14.2|     24.7|

Table: Summary statistics of point-by-point Similarity--Krige surface heat flow estimate differences (mW m$^{-2}$) within each SubMap transect-set domain. {#tbl:est-diff-summary}

\clearpage

|SubMap Set | Region|    n|    Min|   Max| Median|  IQR|  Mean| $\sigma$|
|:----------|------:|----:|------:|-----:|------:|----:|-----:|--------:|
|NPA_SET1   |    NPA|  714|  -58.2|  36.5|  -37.9| 17.1| -32.9|     15.1|
|NPA_SET2   |    NPA|  912|  -38.2|  61.4|  -11.1| 16.5|  -8.3|     15.1|
|NPA_SET3   |    NPA|  758|  -52.2|  39.9|  -27.6| 22.1| -25.8|     15.7|
|NPA_SET4   |    NPA|  755|  -54.5|  33.2|  -32.7| 23.5| -28.8|     15.1|
|NPA_SET5   |    NPA|  912|  -56.1|  44.0|  -24.6| 18.9| -22.3|     15.4|
|NPA_SET6   |    NPA|  401|  -38.6|  39.6|  -10.7| 20.2| -11.8|     14.0|
|NPA_SET7   |    NPA|  650|  -37.2|  50.1|   -3.7| 14.6|  -2.9|     11.3|
|NPA_SET8   |    NPA|  997|  -65.4| 277.8|  -35.4| 24.8| -29.0|     35.7|
|SAM_SET1   |    SAM|  641|  -53.2| 337.7|   23.4| 26.8|  25.4|     29.3|
|SAM_SET2   |    SAM|  981|  -52.8| 111.7|  -21.4| 16.4| -20.0|     18.5|
|SAM_SET3   |    SAM|  661|  -31.3|  46.4|  -10.2| 16.7|  -7.0|     13.2|
|SAM_SET4   |    SAM|  993|  -36.9|  48.1|  -25.1| 14.5| -20.5|     13.3|
|SAM_SET5   |    SAM|  797| -104.0| 117.5|  -26.8| 26.5| -27.6|     23.8|
|SAM_SET6   |    SAM|  348|  -66.8|  33.5|  -27.7| 22.1| -26.8|     17.1|
|SAM_SET7   |    SAM|  584|  -78.9|  36.4|  -29.0| 28.8| -28.4|     21.7|
|SAM_SET8   |    SAM|  631|  -62.9|  61.5|  -31.2| 37.1| -22.7|     25.4|
|SAM_SET9   |    SAM|  328|  -59.3| 142.8|   -9.6| 23.5|  -5.2|     25.6|
|SAM_SET10  |    SAM|  329|  -76.2| 111.7|  -42.8| 18.3| -39.3|     22.1|
|SEA_SET1   |    SEA|  734|  -68.4|  33.7|  -34.5| 17.5| -34.4|     14.0|
|SEA_SET2   |    SEA|  528|  -26.4|  58.8|   -2.0| 14.8|   0.0|     13.4|
|SEA_SET3   |    SEA|  442|  -29.6|  30.4|   -9.8| 14.0|  -8.0|     11.4|
|SEA_SET4   |    SEA|  648|  -33.9|  43.4|   -5.6| 15.4|  -3.6|     12.8|
|SEA_SET5   |    SEA|  502|  -42.9| 221.4|  -10.9| 22.5|  -7.9|     23.1|
|SEA_SET6   |    SEA| 1600|  -35.7| 376.4|    1.0| 21.8|   2.3|     24.2|
|SEA_SET7   |    SEA|  600|  -55.6|  36.7|  -30.2| 18.0| -28.1|     14.5|
|SEA_SET8   |    SEA|  858|  -71.0|  14.7|  -45.0| 19.7| -43.9|     14.2|
|SEA_SET9   |    SEA|  872|  -57.9|  32.7|  -40.1| 21.3| -35.9|     15.3|
|SEA_SET10  |    SEA|  767|  -58.3|  12.4|  -41.2| 13.7| -37.6|     11.5|
|SEA_SET11  |    SEA|  395|  -77.8| 190.7|  -48.8| 27.7| -43.1|     28.8|
|SEA_SET12  |    SEA|  242|  -42.7| 365.5|  -14.4| 19.1| -10.3|     32.5|
|SWP_SET1   |    SWP|  593|  -38.3|  69.4|   -8.8| 24.5|  -7.2|     17.6|
|SWP_SET2   |    SWP|  551|  -80.1| 442.5|  -42.0| 30.1| -27.9|     54.3|
|SWP_SET3   |    SWP|  572|  -50.6| 395.4|  -24.1| 29.2| -16.6|     39.1|
|SWP_SET4   |    SWP|  842|  -68.6| 161.7|  -32.3| 23.2| -32.9|     17.4|

Table: Summary statistics of point-by-point Similarity--Krige estimate-uncertainty differences (mW m$^{-2}$) within each SubMap transect-set domain. {#tbl:sigma-diff-summary}

\clearpage

|SubMap Set | Region|   n|    Min|   Max| Median|   IQR|  Mean| $\sigma$| RMSE|
|:----------|------:|---:|------:|-----:|------:|-----:|-----:|--------:|----:|
|NPA_SET1   |    NPA| 173|  -94.5| 136.6|    5.0|  25.6|   5.3|     29.8| 30.2|
|NPA_SET2   |    NPA| 118|  -63.7| 102.4|    3.3|   7.5|   4.0|     23.7| 23.9|
|NPA_SET3   |    NPA|  48|  -97.1|  42.9|    3.2|  28.8|  -7.4|     32.9| 33.4|
|NPA_SET4   |    NPA|  29|  -97.1|  42.9|    8.6|  23.7|  -1.0|     40.5| 39.8|
|NPA_SET5   |    NPA|  14|  -33.5|  54.5|   25.0|  34.5|  18.8|     28.2| 33.0|
|NPA_SET6   |    NPA|   8|  -33.5|  68.2|    6.6|  27.3|   8.9|     35.2| 34.1|
|NPA_SET7   |    NPA|  17|   -0.1|  68.2|    4.0|  14.3|  15.0|     20.6| 25.0|
|NPA_SET8   |    NPA| 720| -127.3| 122.7|   20.9|  48.7|  19.2|     47.0| 50.7|
|SAM_SET1   |    SAM|  80|  -29.6| 122.9|    3.0|  34.5|  22.4|     36.4| 42.5|
|SAM_SET2   |    SAM|  99| -123.2|  90.9|   11.7|  13.4|  11.7|     29.5| 31.6|
|SAM_SET3   |    SAM|  10|  -33.8|  45.3|    1.9|  12.6|   3.9|     20.8| 20.1|
|SAM_SET4   |    SAM|  31|  -54.2|  85.8|    7.8|  17.2|   2.8|     28.0| 27.7|
|SAM_SET5   |    SAM| 409|  -53.1|  91.0|  -30.8|   0.0| -25.2|     17.5| 30.7|
|SAM_SET6   |    SAM|   2|   -4.5|   9.4|    2.4|   6.9|   2.4|      9.8|  7.4|
|SAM_SET7   |    SAM|  19|  -18.8|  60.7|    0.0|  16.2|   6.7|     25.0| 25.2|
|SAM_SET8   |    SAM|  31|  -55.0|  58.6|    8.5|  16.9|   4.4|     19.2| 19.4|
|SAM_SET9   |    SAM|  38|  -86.1|  72.0|    8.5|  29.2|   4.7|     25.9| 26.0|
|SAM_SET10  |    SAM|  12|  -86.1|  72.0|   -4.5|  25.0|  -1.0|     36.7| 35.1|
|SEA_SET1   |    SEA|  39|  -30.6|  62.0|   -0.9|  12.6|  -0.8|     12.8| 12.7|
|SEA_SET2   |    SEA|  59|  -48.6|  62.0|   -2.3|  16.5|   0.9|     21.4| 21.2|
|SEA_SET3   |    SEA|  13|  -18.4|  59.6|    0.7|   7.3|   2.3|     19.1| 18.5|
|SEA_SET4   |    SEA|  10|  -10.5|  19.3|    5.2|   9.3|   4.5|      9.3|  9.9|
|SEA_SET5   |    SEA|  16|  -84.2|  24.6|    8.5|  14.7|   4.5|     25.3| 24.9|
|SEA_SET6   |    SEA|  41|  -84.2|  75.8|   17.6|  38.5|  22.2|     27.7| 35.2|
|SEA_SET7   |    SEA|  55|  -95.3| 102.8|    1.4|   8.7|  -4.0|     43.4| 43.2|
|SEA_SET8   |    SEA|  63|  -95.3| 132.4|    0.2|  12.8|  -1.2|     49.1| 48.7|
|SEA_SET9   |    SEA| 211|  -94.5| 136.6|    5.0|  23.3|   6.7|     30.0| 30.7|
|SEA_SET10  |    SEA|  69|  -76.6|  87.1|   20.1|  32.1|  10.3|     30.9| 32.4|
|SEA_SET11  |    SEA|  12| -115.2| 123.0|   28.1| 121.2|  31.9|     77.2| 80.5|
|SEA_SET12  |    SEA|   6|  -18.0|  53.0|   -5.8|  33.7|   6.8|     28.2| 26.6|
|SWP_SET1   |    SWP|  21|  -30.8|  27.4|    2.0|  24.9|   6.5|     18.5| 19.2|
|SWP_SET2   |    SWP|  13|  -36.9| 110.4|    2.0|  65.5|  25.2|     43.8| 49.0|
|SWP_SET3   |    SWP|  21|  -12.6| 135.3|   -0.2|  46.6|  24.8|     47.2| 52.3|
|SWP_SET4   |    SWP|  53|  -31.4| 214.8|   -6.8|  44.4|  18.8|     79.7| 81.1|

Table: Summary of Similarity--Observation residuals. {#tbl:sim-residual-summary}

\clearpage

|SubMap Set | Region|   n|    Min|  Max| Median|  IQR|  Mean| $\sigma$| RMSE|
|:----------|------:|---:|------:|----:|------:|----:|-----:|--------:|----:|
|NPA_SET1   |    NPA| 173|  -47.5| 47.2|    2.5|  9.2|   2.2|     10.8| 11.0|
|NPA_SET2   |    NPA| 118|  -17.7| 27.7|   -0.7|  8.0|  -1.7|      5.9|  6.1|
|NPA_SET3   |    NPA|  48|  -40.1| 33.7|    0.0| 11.1|   0.3|     13.2| 13.1|
|NPA_SET4   |    NPA|  29|  -55.4| 75.5|    0.0| 23.8|   2.5|     33.3| 32.8|
|NPA_SET5   |    NPA|  14|  -11.6|  8.8|    0.4|  3.7|   0.4|      5.2|  5.0|
|NPA_SET6   |    NPA|   8|  -18.4|  4.2|    0.2|  6.9|  -3.3|      9.4|  9.4|
|NPA_SET7   |    NPA|  17|   -2.5| 14.1|    0.0|  1.0|   1.3|      3.9|  4.0|
|NPA_SET8   |    NPA| 720| -116.1| 52.8|    2.8| 20.2|   3.0|     22.4| 22.5|
|SAM_SET1   |    SAM|  15|  -52.1| 23.7|    7.5| 35.4|  -7.8|     23.3| 23.8|
|SAM_SET2   |    SAM|  99|  -73.5|  9.5|   -2.4|  3.7|  -6.1|     15.0| 16.1|
|SAM_SET3   |    SAM|  10|   -3.6| 31.7|   -0.3|  1.7|   2.6|     10.3| 10.1|
|SAM_SET4   |    SAM|  16|  -22.4| 82.4|   -1.3|  8.8|   2.8|     23.0| 22.4|
|SAM_SET5   |    SAM| 406|  -21.3|  8.5|  -21.3|  0.0| -18.6|      8.5| 20.4|
|SAM_SET6   |    SAM|   2|   -2.1|  1.6|   -0.2|  1.8|  -0.2|      2.6|  1.9|
|SAM_SET7   |    SAM|  19|  -12.4| 12.7|    0.9|  4.4|  -0.4|      5.0|  4.9|
|SAM_SET8   |    SAM|   5|  -10.1| -0.0|   -1.0|  0.0|  -2.6|      4.2|  4.6|
|SAM_SET9   |    SAM|  38|  -31.0| 12.9|   -0.8|  2.1|  -0.8|      6.5|  6.5|
|SAM_SET10  |    SAM|  12|  -18.0| 21.8|    3.8|  3.8|   2.7|      8.9|  8.9|
|SEA_SET1   |    SEA|  39|  -12.1|  9.1|    0.7|  0.3|   0.1|      3.9|  3.8|
|SEA_SET2   |    SEA|  59|  -23.9| 22.8|   -0.1|  8.3|  -0.3|      7.6|  7.6|
|SEA_SET3   |    SEA|  13|  -11.3|  4.7|   -0.6|  2.3|  -1.2|      4.2|  4.2|
|SEA_SET4   |    SEA|  10|   -1.9|  4.2|   -0.2|  2.5|   0.4|      2.1|  2.0|
|SEA_SET5   |    SEA|  16|  -55.5|  8.1|    4.5|  5.8|   0.4|     15.3| 14.8|
|SEA_SET6   |    SEA|  41|  -55.9|  8.0|    0.7|  6.7|   0.7|      9.6|  9.6|
|SEA_SET7   |    SEA|  55|  -32.7| 24.3|    0.1|  2.3|   0.6|      6.9|  6.9|
|SEA_SET8   |    SEA|  63|  -33.3| 23.4|   -0.0|  4.7|  -0.7|      7.3|  7.3|
|SEA_SET9   |    SEA| 211|  -32.8| 38.8|    0.1|  3.4|   1.6|     10.8| 10.9|
|SEA_SET10  |    SEA|  68|  -21.6| 34.4|    7.8| 15.8|   8.3|     13.2| 15.5|
|SEA_SET11  |    SEA|  12| -102.2| 18.9|  -19.8| 49.7| -23.1|     35.7| 41.3|
|SEA_SET12  |    SEA|   6|   -0.1|  5.2|    1.9|  3.3|   2.1|      2.3|  2.9|
|SWP_SET1   |    SWP|  21|   -8.3|  6.5|    0.0|  4.6|  -1.0|      4.3|  4.3|
|SWP_SET2   |    SWP|  13|  -22.0| 10.8|    0.0|  4.1|  -1.3|      9.8|  9.5|
|SWP_SET3   |    SWP|  13|  -10.7|  7.5|    0.0|  3.8|  -0.9|      5.4|  5.3|
|SWP_SET4   |    SWP|  53|  -27.9| 73.7|   -3.3| 13.5|   7.8|     28.9| 29.7|

Table: Summary of Krige--Observation residuals. {#tbl:krg-residual-summary}

\clearpage

|SubMap Set | Region| Transect Pair| Obs CCF| Sim CCF| Krg CCF| Interpretation|         UID|
|:----------|------:|-------------:|-------:|-------:|-------:|--------------:|-----------:|
|NPA_SET1   |    NPA|        01--02|  0.9010|  0.8348|  0.9309|           cont|  NPA_SET1_1|
|NPA_SET1   |    NPA|        02--03|  0.8690|  0.6605|  0.9063|           cont|  NPA_SET1_2|
|NPA_SET1   |    NPA|        03--04|  0.9669|  0.8551|  0.9915|           cont|  NPA_SET1_3|
|NPA_SET2   |    NPA|        05--06|  0.7716|  0.7888|  0.8284|      ambiguous|  NPA_SET2_1|
|NPA_SET2   |    NPA|        06--07|  0.9440|  0.4985|  0.7257|      ambiguous|  NPA_SET2_2|
|NPA_SET2   |    NPA|        07--08|  0.6004|  0.6623|  0.7351|      ambiguous|  NPA_SET2_3|
|NPA_SET2   |    NPA|        08--09|  0.8244|  0.9111|  0.5681|      ambiguous|  NPA_SET2_4|
|NPA_SET2   |    NPA|        09--10|  0.9992|  0.7348|  0.8282|      ambiguous|  NPA_SET2_5|
|NPA_SET3   |    NPA|        11--12|  0.6212|  0.8454|  0.9004|        discont|  NPA_SET3_1|
|NPA_SET3   |    NPA|        12--13| -0.3759|  0.8164| -0.5146|        discont|  NPA_SET3_2|
|NPA_SET3   |    NPA|        13--14|  0.6720|  0.8880|  0.3961|        discont|  NPA_SET3_3|
|NPA_SET4   |    NPA|        15--16|  0.8410|  0.9433|  0.8622|      ambiguous|  NPA_SET4_1|
|NPA_SET4   |    NPA|        16--17|  0.8848|  0.8056|  0.1941|      ambiguous|  NPA_SET4_2|
|NPA_SET4   |    NPA|        17--18| -0.6787|  0.8482| -0.2082|      ambiguous|  NPA_SET4_3|
|NPA_SET4   |    NPA|        18--19|        |  0.8137|  0.8209|      ambiguous|  NPA_SET4_4|
|NPA_SET4   |    NPA|        19--20|        |  0.8951|  0.9256|      ambiguous|  NPA_SET4_5|
|NPA_SET4   |    NPA|        20--21|        |  0.4706| -0.0442|      ambiguous|  NPA_SET4_6|
|NPA_SET4   |    NPA|        21--22|        |  0.5566|  0.9313|      ambiguous|  NPA_SET4_7|
|NPA_SET5   |    NPA|        23--24|        |  0.9331|  0.7695|      ambiguous|  NPA_SET5_1|
|NPA_SET5   |    NPA|        24--25|        |  0.9712|  0.8871|      ambiguous|  NPA_SET5_2|
|NPA_SET5   |    NPA|        25--26|        |  0.8832|  0.9322|      ambiguous|  NPA_SET5_3|
|NPA_SET5   |    NPA|        26--27|        |  0.9297|  0.8656|      ambiguous|  NPA_SET5_4|
|NPA_SET5   |    NPA|        27--28|        |  0.9617|  0.8968|      ambiguous|  NPA_SET5_5|
|NPA_SET5   |    NPA|        28--29|        |  0.4033|  0.5354|      ambiguous|  NPA_SET5_6|
|NPA_SET5   |    NPA|        29--30|        |  0.9475|  0.1177|      ambiguous|  NPA_SET5_7|
|NPA_SET6   |    NPA|        31--32|        |  0.7617|  0.6158|      ambiguous|  NPA_SET6_1|
|NPA_SET6   |    NPA|        32--33|        |  0.8293|  0.4956|      ambiguous|  NPA_SET6_2|
|NPA_SET6   |    NPA|        33--34|        |  0.8081|  0.9863|      ambiguous|  NPA_SET6_3|
|NPA_SET6   |    NPA|        34--35|        |  0.8751|  0.7702|      ambiguous|  NPA_SET6_4|
|NPA_SET6   |    NPA|        35--36|        |  0.7379| -0.2186|      ambiguous|  NPA_SET6_5|
|NPA_SET7   |    NPA|        37--38|        |  0.6880|  0.9741|      ambiguous|  NPA_SET7_1|
|NPA_SET7   |    NPA|        38--39|        |  0.8933|  0.6305|      ambiguous|  NPA_SET7_2|
|NPA_SET7   |    NPA|        39--40|        |  0.8933|  0.7707|      ambiguous|  NPA_SET7_3|
|NPA_SET7   |    NPA|        40--41|  0.3851|  0.9751|  0.8480|      ambiguous|  NPA_SET7_4|
|NPA_SET7   |    NPA|        41--42|  0.2397|  0.9371|  0.8889|      ambiguous|  NPA_SET7_5|
|NPA_SET8   |    NPA|        43--44|  0.7618|  0.7915|  0.8181|           cont|  NPA_SET8_1|
|NPA_SET8   |    NPA|        44--45|  0.6822|  0.9910|  0.8893|           cont|  NPA_SET8_2|
|NPA_SET8   |    NPA|        45--46|  0.7241|  0.9540|  0.7900|           cont|  NPA_SET8_3|
|NPA_SET8   |    NPA|        46--47|  0.8711|  0.9829|  0.8639|           cont|  NPA_SET8_4|
|SAM_SET1   |    SAM|        01--02|  0.2347|  0.9468|        |      ambiguous|  SAM_SET1_1|
|SAM_SET1   |    SAM|        02--03|  0.7115|  0.9023|        |      ambiguous|  SAM_SET1_2|
|SAM_SET1   |    SAM|        03--04|  0.4599|  0.9686|        |      ambiguous|  SAM_SET1_3|
|SAM_SET1   |    SAM|        04--05|  0.6048|  0.9640|        |      ambiguous|  SAM_SET1_4|
|SAM_SET1   |    SAM|        05--06|  0.9957|  0.8920|        |      ambiguous|  SAM_SET1_5|
|SAM_SET1   |    SAM|        06--07|  0.7077|  0.8148|        |      ambiguous|  SAM_SET1_6|
|SAM_SET1   |    SAM|        07--08|        |  0.9289| -0.3898|      ambiguous|  SAM_SET1_7|
|SAM_SET2   |    SAM|        09--10|  0.9996|  0.3168|  0.9945|        discont|  SAM_SET2_1|
|SAM_SET2   |    SAM|        10--11| -0.2182|  0.9028|  0.2725|        discont|  SAM_SET2_2|
|SAM_SET2   |    SAM|        11--12|  0.7309|  0.9943|  0.8515|        discont|  SAM_SET2_3|
|SAM_SET2   |    SAM|        12--13|        | -0.4955| -0.5377|        discont|  SAM_SET2_4|
|SAM_SET2   |    SAM|        13--14|        |  0.9278|  0.9870|        discont|  SAM_SET2_5|
|SAM_SET2   |    SAM|        14--15|  0.7295|  0.8256|  0.9850|        discont|  SAM_SET2_6|
|SAM_SET3   |    SAM|        16--17|        |  0.9158|  0.7717|      ambiguous|  SAM_SET3_1|
|SAM_SET3   |    SAM|        17--18|        |  0.9164|  0.7233|      ambiguous|  SAM_SET3_2|
|SAM_SET3   |    SAM|        18--19|        |  0.9164|  0.8600|      ambiguous|  SAM_SET3_3|
|SAM_SET3   |    SAM|        19--20|        |  0.7005|  0.8163|      ambiguous|  SAM_SET3_4|
|SAM_SET3   |    SAM|        20--21|        |  0.8965|  0.9237|      ambiguous|  SAM_SET3_5|
|SAM_SET4   |    SAM|        22--23|        |  0.9403| -0.4135|       disagree|  SAM_SET4_1|
|SAM_SET4   |    SAM|        23--24|  0.2381|  0.6678|  0.9640|       disagree|  SAM_SET4_2|
|SAM_SET4   |    SAM|        24--25|  0.7145|  0.4175|  0.1452|       disagree|  SAM_SET4_3|
|SAM_SET4   |    SAM|        25--26|        |  0.9017|  0.0549|       disagree|  SAM_SET4_4|
|SAM_SET4   |    SAM|        26--27|        |  0.0494|  0.9100|       disagree|  SAM_SET4_5|
|SAM_SET4   |    SAM|        27--28|  0.9997| -0.6092|  0.0918|       disagree|  SAM_SET4_6|
|SAM_SET4   |    SAM|        28--29| -0.7331|  0.4914|  0.9671|       disagree|  SAM_SET4_7|
|SAM_SET4   |    SAM|        29--30| -0.7245|  0.7398|  0.9535|       disagree|  SAM_SET4_8|
|SAM_SET4   |    SAM|        30--31| -0.0949|  0.9861|  0.0886|       disagree|  SAM_SET4_9|
|SAM_SET4   |    SAM|        31--32|  0.3099|  0.9125|  0.9246|       disagree| SAM_SET4_10|
|SAM_SET4   |    SAM|        32--33|  0.7557|  0.7240|  0.9905|       disagree| SAM_SET4_11|
|SAM_SET5   |    SAM|        34--35|        |  0.8032|  0.5895|      ambiguous|  SAM_SET5_1|
|SAM_SET5   |    SAM|        35--36|        |  0.9275|  0.7535|      ambiguous|  SAM_SET5_2|
|SAM_SET5   |    SAM|        36--37| -0.4865|  0.9633|  0.7433|      ambiguous|  SAM_SET5_3|
|SAM_SET5   |    SAM|        37--38|        |  0.8464|  0.2030|      ambiguous|  SAM_SET5_4|
|SAM_SET5   |    SAM|        38--39|        |  0.9847|  0.8823|      ambiguous|  SAM_SET5_5|
|SAM_SET5   |    SAM|        39--40|        |  0.8930|  0.7918|      ambiguous|  SAM_SET5_6|
|SAM_SET6   |    SAM|        41--42|        |  0.9270|  0.6326|      ambiguous|  SAM_SET6_1|
|SAM_SET6   |    SAM|        42--43|        |  0.9651|  0.7007|      ambiguous|  SAM_SET6_2|
|SAM_SET6   |    SAM|        43--44|        |  0.1330| -0.1156|      ambiguous|  SAM_SET6_3|
|SAM_SET6   |    SAM|        44--45|        |  0.8651|  0.9116|      ambiguous|  SAM_SET6_4|
|SAM_SET6   |    SAM|        45--46|        |  0.8542|  0.2171|      ambiguous|  SAM_SET6_5|
|SAM_SET7   |    SAM|        47--48|        |  0.8398|  0.7943|      ambiguous|  SAM_SET7_1|
|SAM_SET7   |    SAM|        48--49|        |  0.6972|  0.8143|      ambiguous|  SAM_SET7_2|
|SAM_SET7   |    SAM|        49--50|        |  0.9931|  0.3442|      ambiguous|  SAM_SET7_3|
|SAM_SET7   |    SAM|        50--51|        |  0.7156|  0.9500|      ambiguous|  SAM_SET7_4|
|SAM_SET8   |    SAM|        52--53|        |  0.8651| -0.4936|      ambiguous|  SAM_SET8_1|
|SAM_SET8   |    SAM|        53--54|        |  0.6500|  0.0439|      ambiguous|  SAM_SET8_2|
|SAM_SET8   |    SAM|        54--55|        |  0.5187| -0.4686|      ambiguous|  SAM_SET8_3|
|SAM_SET8   |    SAM|        55--56|        |  0.2707|        |      ambiguous|  SAM_SET8_4|
|SAM_SET9   |    SAM|        57--58|        |  0.4394|  0.5414|      ambiguous|  SAM_SET9_1|
|SAM_SET9   |    SAM|        58--59|        |  0.4475|  0.9802|      ambiguous|  SAM_SET9_2|
|SAM_SET9   |    SAM|        59--60|        |  0.8990|  0.8010|      ambiguous|  SAM_SET9_3|
|SAM_SET9   |    SAM|        60--61|  0.0203|  0.9406| -0.5976|      ambiguous|  SAM_SET9_4|
|SAM_SET10  |    SAM|        62--63|        |  0.9534| -0.1591|      ambiguous| SAM_SET10_1|
|SAM_SET10  |    SAM|        63--64|        |  0.9174|        |      ambiguous| SAM_SET10_2|
|SAM_SET10  |    SAM|        64--65|        |  0.9817|        |      ambiguous| SAM_SET10_3|
|SAM_SET10  |    SAM|        65--66|        |  0.9336|  0.7548|      ambiguous| SAM_SET10_4|
|SEA_SET1   |    SEA|        01--02|        |  0.9817|  0.9126|      ambiguous|  SEA_SET1_1|
|SEA_SET1   |    SEA|        02--03|        |  0.9487|  0.8523|      ambiguous|  SEA_SET1_2|
|SEA_SET1   |    SEA|        03--04|        |  0.9949|  0.9259|      ambiguous|  SEA_SET1_3|
|SEA_SET1   |    SEA|        04--05|        |  0.8877|  0.7198|      ambiguous|  SEA_SET1_4|
|SEA_SET1   |    SEA|        05--06|        |  0.9504|  0.8020|      ambiguous|  SEA_SET1_5|
|SEA_SET1   |    SEA|        06--07|        |  1.0000|  0.8429|      ambiguous|  SEA_SET1_6|
|SEA_SET2   |    SEA|        08--09|  0.8712|  0.9376|  0.9399|           cont|  SEA_SET2_1|
|SEA_SET2   |    SEA|        09--10|  0.9998|  0.8377|  0.8782|           cont|  SEA_SET2_2|
|SEA_SET2   |    SEA|        10--11|  0.5615|  0.9260|  0.7037|           cont|  SEA_SET2_3|
|SEA_SET2   |    SEA|        11--12|  0.7140|  0.9686|  0.8817|           cont|  SEA_SET2_4|
|SEA_SET3   |    SEA|        13--14|  0.5494|  0.7594| -0.0050|      ambiguous|  SEA_SET3_1|
|SEA_SET3   |    SEA|        14--15|        |  0.7639|  0.4423|      ambiguous|  SEA_SET3_2|
|SEA_SET3   |    SEA|        15--16|        |  0.7457|  0.4413|      ambiguous|  SEA_SET3_3|
|SEA_SET3   |    SEA|        16--17|        |  0.7307|  0.8026|      ambiguous|  SEA_SET3_4|
|SEA_SET4   |    SEA|        18--19|        |  0.9719|  0.8292|      ambiguous|  SEA_SET4_1|
|SEA_SET4   |    SEA|        19--20|        |  0.4772|  0.0575|      ambiguous|  SEA_SET4_2|
|SEA_SET4   |    SEA|        20--21|        |  0.7321|  0.8876|      ambiguous|  SEA_SET4_3|
|SEA_SET4   |    SEA|        21--22|        |  0.9134|  0.9386|      ambiguous|  SEA_SET4_4|
|SEA_SET4   |    SEA|        22--29|        | -0.5155|  0.0563|      ambiguous|  SEA_SET4_5|
|SEA_SET4   |    SEA|        29--30|        |  0.9222|  0.0918|      ambiguous|  SEA_SET4_6|
|SEA_SET5   |    SEA|        23--24|        |  0.8047|  0.2620|      ambiguous|  SEA_SET5_1|
|SEA_SET5   |    SEA|        24--25|        |  0.8794|  0.6410|      ambiguous|  SEA_SET5_2|
|SEA_SET5   |    SEA|        25--26|        |  0.9660|  0.9623|      ambiguous|  SEA_SET5_3|
|SEA_SET5   |    SEA|        26--27|  0.7611|  0.8518|  0.7433|      ambiguous|  SEA_SET5_4|
|SEA_SET5   |    SEA|        27--28|  0.8694|  0.9006|  0.4252|      ambiguous|  SEA_SET5_5|
|SEA_SET5   |    SEA|        28--31|        | -0.6596| -0.0836|      ambiguous|  SEA_SET5_6|
|SEA_SET5   |    SEA|        31--32|        | -0.2257| -0.0601|      ambiguous|  SEA_SET5_7|
|SEA_SET5   |    SEA|        32--33|  0.7571|  0.9322|  0.7295|      ambiguous|  SEA_SET5_8|
|SEA_SET5   |    SEA|        33--34|        |  0.9583| -0.3179|      ambiguous|  SEA_SET5_9|
|SEA_SET6   |    SEA|        35--36|        |  0.9889|  0.4148|      ambiguous|  SEA_SET6_1|
|SEA_SET6   |    SEA|        36--37|        |  0.4997|  0.2648|      ambiguous|  SEA_SET6_2|
|SEA_SET6   |    SEA|        37--38| -0.1408|  0.5026| -0.4231|      ambiguous|  SEA_SET6_3|
|SEA_SET6   |    SEA|        38--39|  0.9306|  0.8325| -0.3113|      ambiguous|  SEA_SET6_4|
|SEA_SET6   |    SEA|        39--40|        | -0.7013|  0.6490|      ambiguous|  SEA_SET6_5|
|SEA_SET6   |    SEA|        40--41|        |  0.9999|  0.4365|      ambiguous|  SEA_SET6_6|
|SEA_SET6   |    SEA|        41--42|        |  0.8578|  0.2249|      ambiguous|  SEA_SET6_7|
|SEA_SET6   |    SEA|        42--43|        |  0.8588|  0.5046|      ambiguous|  SEA_SET6_8|
|SEA_SET6   |    SEA|        43--44|        |  1.0000|  0.9331|      ambiguous|  SEA_SET6_9|
|SEA_SET6   |    SEA|        44--45|        |  1.0000|  0.8779|      ambiguous| SEA_SET6_10|
|SEA_SET6   |    SEA|        45--46|        |  0.9419|  0.3738|      ambiguous| SEA_SET6_11|
|SEA_SET6   |    SEA|        46--47|        |  0.8513| -0.1997|      ambiguous| SEA_SET6_12|
|SEA_SET6   |    SEA|        47--48|        |  0.8491|  0.6518|      ambiguous| SEA_SET6_13|
|SEA_SET6   |    SEA|        48--49|        |  0.5115|  0.2840|      ambiguous| SEA_SET6_14|
|SEA_SET6   |    SEA|        49--50|        |  0.1779|  0.7533|      ambiguous| SEA_SET6_15|
|SEA_SET6   |    SEA|        50--51|        |  0.4431|  0.2947|      ambiguous| SEA_SET6_16|
|SEA_SET6   |    SEA|        51--52|        |  0.8946| -0.5184|      ambiguous| SEA_SET6_17|
|SEA_SET6   |    SEA|        52--53|        |  0.9103|  0.5435|      ambiguous| SEA_SET6_18|
|SEA_SET6   |    SEA|        53--54|        |  0.8629|  0.4541|      ambiguous| SEA_SET6_19|
|SEA_SET7   |    SEA|        55--56|        |  0.7875|  0.6749|        discont|  SEA_SET7_1|
|SEA_SET7   |    SEA|        56--57|  0.6740|  0.5647|  0.7123|        discont|  SEA_SET7_2|
|SEA_SET7   |    SEA|        57--58|  0.3621|  0.5252|  0.2147|        discont|  SEA_SET7_3|
|SEA_SET7   |    SEA|        58--59|  0.3122| -0.2729|  0.8018|        discont|  SEA_SET7_4|
|SEA_SET7   |    SEA|        59--60| -0.2346| -0.3734|  0.7516|        discont|  SEA_SET7_5|
|SEA_SET8   |    SEA|        61--62|  0.9703|  0.7665|  0.3518|      ambiguous|  SEA_SET8_1|
|SEA_SET8   |    SEA|        62--63|  0.4639|  0.3858|  0.5928|      ambiguous|  SEA_SET8_2|
|SEA_SET8   |    SEA|        63--64|  0.8009|  0.9333|  0.4346|      ambiguous|  SEA_SET8_3|
|SEA_SET8   |    SEA|        64--65|  0.9691|  0.9076|  0.9266|      ambiguous|  SEA_SET8_4|
|SEA_SET8   |    SEA|        65--66|  0.4269|  0.9736|  0.8255|      ambiguous|  SEA_SET8_5|
|SEA_SET9   |    SEA|        67--68|  0.2741| -0.4037|  0.3003|        discont|  SEA_SET9_1|
|SEA_SET9   |    SEA|        68--69| -0.3134|  0.5301| -0.1464|        discont|  SEA_SET9_2|
|SEA_SET9   |    SEA|        69--70| -0.3025|  0.7943| -0.1980|        discont|  SEA_SET9_3|
|SEA_SET9   |    SEA|        70--71|  0.9998|  0.1303|  0.9999|        discont|  SEA_SET9_4|
|SEA_SET9   |    SEA|        71--72|  0.5031|  0.4402|  0.8693|        discont|  SEA_SET9_5|
|SEA_SET10  |    SEA|        73--74| -0.2328|  0.9282|  0.3381|       disagree| SEA_SET10_1|
|SEA_SET10  |    SEA|        74--75|  0.8831|  0.8154| -0.0027|       disagree| SEA_SET10_2|
|SEA_SET10  |    SEA|        75--76|  0.7915|  0.8866|  0.6436|       disagree| SEA_SET10_3|
|SEA_SET10  |    SEA|        76--77| -0.0757|  0.8605|  0.1541|       disagree| SEA_SET10_4|
|SEA_SET11  |    SEA|        78--79| -0.2744|  0.9739|  0.1411|      ambiguous| SEA_SET11_1|
|SEA_SET11  |    SEA|        79--80|  0.3636|  1.0000|  0.7268|      ambiguous| SEA_SET11_2|
|SEA_SET11  |    SEA|        80--81|  0.4639|  0.9945|  0.8978|      ambiguous| SEA_SET11_3|
|SEA_SET11  |    SEA|        81--82|  0.9918|  0.9274|  0.9278|      ambiguous| SEA_SET11_4|
|SEA_SET11  |    SEA|        82--83|  0.9459|  0.9178|  0.9465|      ambiguous| SEA_SET11_5|
|SEA_SET11  |    SEA|        83--84|        |  0.8756|  0.9284|      ambiguous| SEA_SET11_6|
|SEA_SET11  |    SEA|        84--85|        |  0.9321|  0.7678|      ambiguous| SEA_SET11_7|
|SEA_SET11  |    SEA|        85--86|        |  0.8987|  0.7279|      ambiguous| SEA_SET11_8|
|SEA_SET12  |    SEA|        87--88|        |  0.9442|  0.4201|      ambiguous| SEA_SET12_1|
|SEA_SET12  |    SEA|        88--89|        | -0.0983|  0.8384|      ambiguous| SEA_SET12_2|
|SEA_SET12  |    SEA|        89--90|        |  0.4235|  0.0093|      ambiguous| SEA_SET12_3|
|SWP_SET1   |    SWP|        01--02|        | -0.2606|  0.8677|      ambiguous|  SWP_SET1_1|
|SWP_SET1   |    SWP|        02--03|        |  0.2007|  0.5396|      ambiguous|  SWP_SET1_2|
|SWP_SET1   |    SWP|        03--04|        |  0.3468| -0.4069|      ambiguous|  SWP_SET1_3|
|SWP_SET1   |    SWP|        04--05|        |  0.9849|  0.9233|      ambiguous|  SWP_SET1_4|
|SWP_SET1   |    SWP|        05--06|        |  0.8860|  0.7946|      ambiguous|  SWP_SET1_5|
|SWP_SET1   |    SWP|        06--07|        |  0.8677|  0.8479|      ambiguous|  SWP_SET1_6|
|SWP_SET1   |    SWP|        07--08|        |  0.8748|        |      ambiguous|  SWP_SET1_7|
|SWP_SET1   |    SWP|        08--09|        |  0.8928|        |      ambiguous|  SWP_SET1_8|
|SWP_SET2   |    SWP|        10--11|        |  0.6515| -0.1297|      ambiguous|  SWP_SET2_1|
|SWP_SET2   |    SWP|        11--12|        |  0.9753|  0.7179|      ambiguous|  SWP_SET2_2|
|SWP_SET2   |    SWP|        12--13|  0.6894|  0.8505|  0.1682|      ambiguous|  SWP_SET2_3|
|SWP_SET2   |    SWP|        13--14|        |  0.9749|  0.7456|      ambiguous|  SWP_SET2_4|
|SWP_SET2   |    SWP|        14--15|        |  0.9871|  0.2872|      ambiguous|  SWP_SET2_5|
|SWP_SET2   |    SWP|        15--16|        |  0.8662|  0.9259|      ambiguous|  SWP_SET2_6|
|SWP_SET3   |    SWP|        17--18|        |  0.7457|  0.3938|      ambiguous|  SWP_SET3_1|
|SWP_SET3   |    SWP|        18--19|        |  0.9028|  0.9126|      ambiguous|  SWP_SET3_2|
|SWP_SET3   |    SWP|        19--20|        |  0.9278|  0.8014|      ambiguous|  SWP_SET3_3|
|SWP_SET3   |    SWP|        20--21|        |  0.8312|  0.6610|      ambiguous|  SWP_SET3_4|
|SWP_SET3   |    SWP|        21--22|        |  0.7680| -0.3102|      ambiguous|  SWP_SET3_5|
|SWP_SET3   |    SWP|        22--23|        |  0.8445|  0.0127|      ambiguous|  SWP_SET3_6|
|SWP_SET3   |    SWP|        23--24|        |  0.9485|        |      ambiguous|  SWP_SET3_7|
|SWP_SET3   |    SWP|        24--25|        |  0.8443|        |      ambiguous|  SWP_SET3_8|
|SWP_SET3   |    SWP|        25--26|        |  0.8983|        |      ambiguous|  SWP_SET3_9|
|SWP_SET4   |    SWP|        27--28|        |  0.8630|  0.7999|           cont|  SWP_SET4_1|
|SWP_SET4   |    SWP|        28--29|        |  0.8573|  0.7916|           cont|  SWP_SET4_2|
|SWP_SET4   |    SWP|        29--30|  0.9687|  0.9311|  0.9405|           cont|  SWP_SET4_3|
|SWP_SET4   |    SWP|        30--31|  0.9762|  0.8963|  0.8842|           cont|  SWP_SET4_4|
|SWP_SET4   |    SWP|        31--32|  0.8895|  0.7782|  0.8050|           cont|  SWP_SET4_5|

Table: Global cross-correlations between adjacent SubMap transect pairs by method (observations, Similarity, and Krige). A total of 201 adjacent pairs are listed across 34 transect sets. Empty cells indicate insufficient projected data to compute a correlation. NPA: North Pacific--Aleutians; SAM: South America; SEA: Southeast Asia; SWP: Southwest Pacific. {#tbl:global-ccf-summary}

\clearpage

|SubMap Set | Region| Transect ID| Obs--Sim CCF| Obs--Krg CCF| Sim--Krg CCF| Interpretation|
|:----------|------:|-----------:|------------:|------------:|------------:|--------------:|
|NPA_SET1   |    NPA|       NPA01|       0.9457|       0.9701|       0.9652|           cont|
|NPA_SET1   |    NPA|       NPA02|       0.8459|       0.9645|       0.8332|           cont|
|NPA_SET1   |    NPA|       NPA03|       0.9820|       0.9655|       0.9650|           cont|
|NPA_SET1   |    NPA|       NPA04|       0.9064|       0.9816|       0.9275|           cont|
|NPA_SET2   |    NPA|       NPA05|       0.9137|       0.9876|       0.9308|      ambiguous|
|NPA_SET2   |    NPA|       NPA06|       0.9617|       0.9804|       0.9858|      ambiguous|
|NPA_SET2   |    NPA|       NPA07|       0.7912|       0.9486|       0.8466|      ambiguous|
|NPA_SET2   |    NPA|       NPA08|       0.9831|       0.7543|       0.7369|      ambiguous|
|NPA_SET2   |    NPA|       NPA09|       0.7102|       0.2962|       0.7341|      ambiguous|
|NPA_SET2   |    NPA|       NPA10|       0.3963|       0.6138|       0.8527|      ambiguous|
|NPA_SET3   |    NPA|       NPA11|       0.7670|       0.7829|       0.9447|        discont|
|NPA_SET3   |    NPA|       NPA12|       0.8207|       0.8656|       0.9696|        discont|
|NPA_SET3   |    NPA|       NPA13|      -0.5772|       0.9233|      -0.5924|        discont|
|NPA_SET3   |    NPA|       NPA14|      -0.1029|       0.2792|       0.1893|        discont|
|NPA_SET4   |    NPA|       NPA15|       0.6190|       0.9804|       0.6867|      ambiguous|
|NPA_SET4   |    NPA|       NPA16|       0.1766|       0.6231|       0.5513|      ambiguous|
|NPA_SET4   |    NPA|       NPA17|      -0.1718|       0.1778|       0.8616|      ambiguous|
|NPA_SET4   |    NPA|       NPA18|       0.7694|      -0.1663|      -0.5632|      ambiguous|
|NPA_SET4   |    NPA|       NPA19|             |             |      -0.4364|      ambiguous|
|NPA_SET4   |    NPA|       NPA20|       0.2926|       0.5900|      -0.0528|      ambiguous|
|NPA_SET4   |    NPA|       NPA21|             |             |       0.6479|      ambiguous|
|NPA_SET4   |    NPA|       NPA22|       0.5902|       0.7096|       0.5113|      ambiguous|
|NPA_SET5   |    NPA|       NPA23|             |             |       0.5616|      ambiguous|
|NPA_SET5   |    NPA|       NPA24|       0.6242|       0.9381|       0.7263|      ambiguous|
|NPA_SET5   |    NPA|       NPA25|             |             |       0.8198|      ambiguous|
|NPA_SET5   |    NPA|       NPA26|             |             |       0.6660|      ambiguous|
|NPA_SET5   |    NPA|       NPA27|      -0.0320|       0.3471|       0.6840|      ambiguous|
|NPA_SET5   |    NPA|       NPA28|             |             |       0.8036|      ambiguous|
|NPA_SET5   |    NPA|       NPA29|             |             |      -0.2032|      ambiguous|
|NPA_SET5   |    NPA|       NPA30|             |             |      -0.1569|      ambiguous|
|NPA_SET6   |    NPA|       NPA31|             |             |       0.5554|      ambiguous|
|NPA_SET6   |    NPA|       NPA32|             |             |       0.5042|      ambiguous|
|NPA_SET6   |    NPA|       NPA33|             |             |       0.4619|      ambiguous|
|NPA_SET6   |    NPA|       NPA34|             |             |       0.0882|      ambiguous|
|NPA_SET6   |    NPA|       NPA35|             |             |      -0.2175|      ambiguous|
|NPA_SET6   |    NPA|       NPA36|             |             |       0.4408|      ambiguous|
|NPA_SET7   |    NPA|       NPA37|             |             |      -0.0578|      ambiguous|
|NPA_SET7   |    NPA|       NPA38|             |             |      -0.0110|      ambiguous|
|NPA_SET7   |    NPA|       NPA39|             |             |       0.4031|      ambiguous|
|NPA_SET7   |    NPA|       NPA40|       0.7272|       0.3152|       0.6340|      ambiguous|
|NPA_SET7   |    NPA|       NPA41|       0.7745|       0.9477|       0.6149|      ambiguous|
|NPA_SET7   |    NPA|       NPA42|       0.4793|       0.3326|       0.3373|      ambiguous|
|NPA_SET8   |    NPA|       NPA43|       0.8055|       0.8297|       0.7887|           cont|
|NPA_SET8   |    NPA|       NPA44|       0.6380|       0.7963|       0.8227|           cont|
|NPA_SET8   |    NPA|       NPA45|       0.9536|       0.9870|       0.9356|           cont|
|NPA_SET8   |    NPA|       NPA46|       0.8586|       0.9717|       0.8277|           cont|
|NPA_SET8   |    NPA|       NPA47|       0.9827|       0.9800|       0.9825|           cont|
|SAM_SET1   |    SAM|       SAM01|       0.7564|             |             |      ambiguous|
|SAM_SET1   |    SAM|       SAM02|       0.8463|             |             |      ambiguous|
|SAM_SET1   |    SAM|       SAM03|       0.2119|             |             |      ambiguous|
|SAM_SET1   |    SAM|       SAM04|       0.9731|             |             |      ambiguous|
|SAM_SET1   |    SAM|       SAM05|       0.8357|             |             |      ambiguous|
|SAM_SET1   |    SAM|       SAM06|       0.9008|             |             |      ambiguous|
|SAM_SET1   |    SAM|       SAM07|       0.8454|       0.5182|       0.0716|      ambiguous|
|SAM_SET1   |    SAM|       SAM08|             |             |       0.3872|      ambiguous|
|SAM_SET2   |    SAM|       SAM09|      -0.3027|       0.9516|      -0.2669|        discont|
|SAM_SET2   |    SAM|       SAM10|      -0.3678|       0.9503|      -0.2346|        discont|
|SAM_SET2   |    SAM|       SAM11|       0.8343|       0.6960|       0.4488|        discont|
|SAM_SET2   |    SAM|       SAM12|       0.7072|       0.7581|       0.8815|        discont|
|SAM_SET2   |    SAM|       SAM13|             |             |       0.6317|        discont|
|SAM_SET2   |    SAM|       SAM14|       0.6902|       0.6628|       0.7651|        discont|
|SAM_SET2   |    SAM|       SAM15|       0.8912|       0.9543|       0.8457|        discont|
|SAM_SET3   |    SAM|       SAM16|             |             |       0.7744|      ambiguous|
|SAM_SET3   |    SAM|       SAM17|       0.8920|       0.8829|       0.7397|      ambiguous|
|SAM_SET3   |    SAM|       SAM18|             |             |       0.7005|      ambiguous|
|SAM_SET3   |    SAM|       SAM19|       0.5780|       0.9036|       0.6089|      ambiguous|
|SAM_SET3   |    SAM|       SAM20|             |             |       0.7908|      ambiguous|
|SAM_SET3   |    SAM|       SAM21|       0.5654|       0.4491|       0.7920|      ambiguous|
|SAM_SET4   |    SAM|       SAM22|             |             |      -0.3608|       disagree|
|SAM_SET4   |    SAM|       SAM23|       0.4300|      -0.1459|       0.6296|       disagree|
|SAM_SET4   |    SAM|       SAM24|       0.5525|       0.7429|       0.4194|       disagree|
|SAM_SET4   |    SAM|       SAM25|       0.1070|       0.3173|      -0.0014|       disagree|
|SAM_SET4   |    SAM|       SAM26|             |             |       0.7503|       disagree|
|SAM_SET4   |    SAM|       SAM27|       0.7353|      -0.2532|       0.1094|       disagree|
|SAM_SET4   |    SAM|       SAM28|      -0.6988|      -0.1395|       0.4296|       disagree|
|SAM_SET4   |    SAM|       SAM29|       0.6272|       0.5074|       0.7451|       disagree|
|SAM_SET4   |    SAM|       SAM30|      -0.0818|      -0.0782|       0.9242|       disagree|
|SAM_SET4   |    SAM|       SAM31|       0.7751|       0.1135|       0.1566|       disagree|
|SAM_SET4   |    SAM|       SAM32|       0.6134|      -0.0526|      -0.3557|       disagree|
|SAM_SET4   |    SAM|       SAM33|       0.7610|       0.0793|      -0.1768|       disagree|
|SAM_SET5   |    SAM|       SAM34|             |             |       0.9495|      ambiguous|
|SAM_SET5   |    SAM|       SAM35|             |             |       0.7055|      ambiguous|
|SAM_SET5   |    SAM|       SAM36|       0.3212|       0.8042|       0.1936|      ambiguous|
|SAM_SET5   |    SAM|       SAM37|       0.2712|       0.2851|       0.4608|      ambiguous|
|SAM_SET5   |    SAM|       SAM38|             |             |       0.9759|      ambiguous|
|SAM_SET5   |    SAM|       SAM39|             |             |       0.8585|      ambiguous|
|SAM_SET5   |    SAM|       SAM40|             |             |       0.7202|      ambiguous|
|SAM_SET6   |    SAM|       SAM41|             |             |      -0.1156|      ambiguous|
|SAM_SET6   |    SAM|       SAM42|       0.9173|       0.6389|       0.5922|      ambiguous|
|SAM_SET6   |    SAM|       SAM43|             |             |       0.9146|      ambiguous|
|SAM_SET6   |    SAM|       SAM44|             |             |       0.5839|      ambiguous|
|SAM_SET6   |    SAM|       SAM45|             |             |       0.6068|      ambiguous|
|SAM_SET6   |    SAM|       SAM46|             |             |       0.5459|      ambiguous|
|SAM_SET7   |    SAM|       SAM47|             |             |       0.7575|      ambiguous|
|SAM_SET7   |    SAM|       SAM48|       0.8583|       0.9725|       0.9102|      ambiguous|
|SAM_SET7   |    SAM|       SAM49|             |             |       0.9790|      ambiguous|
|SAM_SET7   |    SAM|       SAM50|             |             |       0.2543|      ambiguous|
|SAM_SET7   |    SAM|       SAM51|             |             |      -0.2034|      ambiguous|
|SAM_SET8   |    SAM|       SAM52|       0.5672|       0.3409|       0.8438|      ambiguous|
|SAM_SET8   |    SAM|       SAM53|             |             |      -0.2799|      ambiguous|
|SAM_SET8   |    SAM|       SAM54|             |             |       0.8593|      ambiguous|
|SAM_SET8   |    SAM|       SAM55|             |             |       0.2744|      ambiguous|
|SAM_SET8   |    SAM|       SAM56|      -0.2678|             |             |      ambiguous|
|SAM_SET9   |    SAM|       SAM57|       0.4756|      -0.0145|       0.2302|      ambiguous|
|SAM_SET9   |    SAM|       SAM58|             |             |       0.5710|      ambiguous|
|SAM_SET9   |    SAM|       SAM59|             |             |       0.6097|      ambiguous|
|SAM_SET9   |    SAM|       SAM60|       0.9287|       0.7246|       0.7917|      ambiguous|
|SAM_SET9   |    SAM|       SAM61|       0.1376|       0.0627|      -0.2678|      ambiguous|
|SAM_SET10  |    SAM|       SAM62|       0.6794|       0.7779|       0.7000|      ambiguous|
|SAM_SET10  |    SAM|       SAM63|             |             |      -0.6913|      ambiguous|
|SAM_SET10  |    SAM|       SAM64|             |             |             |      ambiguous|
|SAM_SET10  |    SAM|       SAM65|             |             |       0.3309|      ambiguous|
|SAM_SET10  |    SAM|       SAM66|             |             |      -0.3219|      ambiguous|
|SEA_SET1   |    SEA|       SEA01|             |             |       0.7304|      ambiguous|
|SEA_SET1   |    SEA|       SEA02|             |             |       0.8599|      ambiguous|
|SEA_SET1   |    SEA|       SEA03|       0.8648|       0.9196|       0.7549|      ambiguous|
|SEA_SET1   |    SEA|       SEA04|             |             |       0.8513|      ambiguous|
|SEA_SET1   |    SEA|       SEA05|             |             |       0.9735|      ambiguous|
|SEA_SET1   |    SEA|       SEA06|             |             |       0.7622|      ambiguous|
|SEA_SET1   |    SEA|       SEA07|       0.8307|       0.9027|       0.9753|      ambiguous|
|SEA_SET2   |    SEA|       SEA08|       0.7543|       0.9673|       0.7228|           cont|
|SEA_SET2   |    SEA|       SEA09|       0.6628|       0.8383|       0.8674|           cont|
|SEA_SET2   |    SEA|       SEA10|       0.8357|       0.6977|       0.8030|           cont|
|SEA_SET2   |    SEA|       SEA11|       0.8812|       0.6779|       0.6360|           cont|
|SEA_SET2   |    SEA|       SEA12|       0.5061|       0.7133|       0.7062|           cont|
|SEA_SET3   |    SEA|       SEA13|      -0.4912|       0.1889|       0.1635|      ambiguous|
|SEA_SET3   |    SEA|       SEA14|       0.0666|       0.8216|       0.4424|      ambiguous|
|SEA_SET3   |    SEA|       SEA15|             |             |       0.1331|      ambiguous|
|SEA_SET3   |    SEA|       SEA16|             |             |       0.7740|      ambiguous|
|SEA_SET3   |    SEA|       SEA17|       0.8551|       0.8553|       1.0000|      ambiguous|
|SEA_SET4   |    SEA|       SEA18|             |             |       0.8786|      ambiguous|
|SEA_SET4   |    SEA|       SEA19|             |             |       0.7504|      ambiguous|
|SEA_SET4   |    SEA|       SEA20|       0.0625|       0.6935|       0.6860|      ambiguous|
|SEA_SET4   |    SEA|       SEA21|             |             |       0.1146|      ambiguous|
|SEA_SET4   |    SEA|       SEA22|             |             |       0.3000|      ambiguous|
|SEA_SET4   |    SEA|       SEA29|             |             |       0.5415|      ambiguous|
|SEA_SET4   |    SEA|       SEA30|             |             |       0.1659|      ambiguous|
|SEA_SET5   |    SEA|       SEA23|             |             |       0.8532|      ambiguous|
|SEA_SET5   |    SEA|       SEA24|             |             |       0.6150|      ambiguous|
|SEA_SET5   |    SEA|       SEA25|             |             |       0.4731|      ambiguous|
|SEA_SET5   |    SEA|       SEA26|       0.1991|       0.9475|       0.2132|      ambiguous|
|SEA_SET5   |    SEA|       SEA27|       0.8030|       0.9774|       0.7216|      ambiguous|
|SEA_SET5   |    SEA|       SEA28|       0.8865|       0.6455|       0.4271|      ambiguous|
|SEA_SET5   |    SEA|       SEA31|             |             |       0.6269|      ambiguous|
|SEA_SET5   |    SEA|       SEA32|       0.1411|       0.6309|       0.4366|      ambiguous|
|SEA_SET5   |    SEA|       SEA33|       0.1335|       0.8378|       0.4005|      ambiguous|
|SEA_SET5   |    SEA|       SEA34|             |             |       0.2727|      ambiguous|
|SEA_SET6   |    SEA|       SEA35|       0.8347|       0.4133|       0.2402|      ambiguous|
|SEA_SET6   |    SEA|       SEA36|             |             |      -0.1360|      ambiguous|
|SEA_SET6   |    SEA|       SEA37|       0.4544|       0.7972|       0.4264|      ambiguous|
|SEA_SET6   |    SEA|       SEA38|       0.9646|       0.7199|       0.7994|      ambiguous|
|SEA_SET6   |    SEA|       SEA39|       0.9505|      -0.3311|      -0.3330|      ambiguous|
|SEA_SET6   |    SEA|       SEA40|             |             |       0.5752|      ambiguous|
|SEA_SET6   |    SEA|       SEA41|             |             |       0.9480|      ambiguous|
|SEA_SET6   |    SEA|       SEA42|             |             |      -0.0445|      ambiguous|
|SEA_SET6   |    SEA|       SEA43|             |             |       0.6421|      ambiguous|
|SEA_SET6   |    SEA|       SEA44|             |             |       0.7564|      ambiguous|
|SEA_SET6   |    SEA|       SEA45|       0.6287|       0.9350|       0.5760|      ambiguous|
|SEA_SET6   |    SEA|       SEA46|             |             |      -0.4600|      ambiguous|
|SEA_SET6   |    SEA|       SEA47|             |             |       0.8363|      ambiguous|
|SEA_SET6   |    SEA|       SEA48|             |             |       0.6440|      ambiguous|
|SEA_SET6   |    SEA|       SEA49|             |             |       0.4007|      ambiguous|
|SEA_SET6   |    SEA|       SEA50|             |             |       0.0719|      ambiguous|
|SEA_SET6   |    SEA|       SEA51|             |             |      -0.3830|      ambiguous|
|SEA_SET6   |    SEA|       SEA52|             |             |       0.8659|      ambiguous|
|SEA_SET6   |    SEA|       SEA53|      -0.0634|       0.0539|       0.6589|      ambiguous|
|SEA_SET6   |    SEA|       SEA54|             |             |       0.7867|      ambiguous|
|SEA_SET7   |    SEA|       SEA55|             |             |       0.8346|        discont|
|SEA_SET7   |    SEA|       SEA56|       0.4920|       0.8516|       0.4490|        discont|
|SEA_SET7   |    SEA|       SEA57|       0.7020|       0.7775|       0.8716|        discont|
|SEA_SET7   |    SEA|       SEA58|       0.3261|       0.9526|       0.2604|        discont|
|SEA_SET7   |    SEA|       SEA59|       0.9947|       0.5096|       0.4784|        discont|
|SEA_SET7   |    SEA|       SEA60|       0.5792|       0.8077|       0.3622|        discont|
|SEA_SET8   |    SEA|       SEA61|       0.8695|       0.6613|       0.7076|      ambiguous|
|SEA_SET8   |    SEA|       SEA62|       0.9974|       0.4187|       0.4183|      ambiguous|
|SEA_SET8   |    SEA|       SEA63|       0.7189|       0.4970|       0.6815|      ambiguous|
|SEA_SET8   |    SEA|       SEA64|       0.8157|       0.9279|       0.8930|      ambiguous|
|SEA_SET8   |    SEA|       SEA65|       0.9261|       0.9408|       0.9300|      ambiguous|
|SEA_SET8   |    SEA|       SEA66|       0.3662|       0.3348|       0.8766|      ambiguous|
|SEA_SET9   |    SEA|       SEA67|       0.0661|       0.9125|      -0.1500|        discont|
|SEA_SET9   |    SEA|       SEA68|       0.6506|       0.9929|       0.6739|        discont|
|SEA_SET9   |    SEA|       SEA69|       0.8730|       0.8859|       0.8205|        discont|
|SEA_SET9   |    SEA|       SEA70|       0.1356|       0.9999|       0.1307|        discont|
|SEA_SET9   |    SEA|       SEA71|       0.9993|       0.9993|       1.0000|        discont|
|SEA_SET9   |    SEA|       SEA72|       0.9761|       0.6954|       0.6365|        discont|
|SEA_SET10  |    SEA|       SEA73|      -0.2396|      -0.0616|       0.9107|       disagree|
|SEA_SET10  |    SEA|       SEA74|       0.8341|       0.3147|       0.5913|       disagree|
|SEA_SET10  |    SEA|       SEA75|       0.7379|       0.9967|       0.7108|       disagree|
|SEA_SET10  |    SEA|       SEA76|       0.8907|       0.9579|       0.8692|       disagree|
|SEA_SET10  |    SEA|       SEA77|       0.1380|       0.8591|       0.2473|       disagree|
|SEA_SET11  |    SEA|       SEA78|      -0.3754|      -0.1769|       0.5397|      ambiguous|
|SEA_SET11  |    SEA|       SEA79|       0.9535|       0.8319|       0.8162|      ambiguous|
|SEA_SET11  |    SEA|       SEA80|       0.2400|       0.7085|       0.5624|      ambiguous|
|SEA_SET11  |    SEA|       SEA81|       0.8616|       0.8698|       0.7050|      ambiguous|
|SEA_SET11  |    SEA|       SEA82|       0.7141|       0.9927|       0.7373|      ambiguous|
|SEA_SET11  |    SEA|       SEA83|       0.9596|       0.9944|       0.9714|      ambiguous|
|SEA_SET11  |    SEA|       SEA84|             |             |       0.8991|      ambiguous|
|SEA_SET11  |    SEA|       SEA85|       0.8291|       0.9771|       0.7532|      ambiguous|
|SEA_SET11  |    SEA|       SEA86|             |             |       0.8139|      ambiguous|
|SEA_SET12  |    SEA|       SEA87|       0.8865|       0.2781|       0.4452|      ambiguous|
|SEA_SET12  |    SEA|       SEA88|             |             |       0.6050|      ambiguous|
|SEA_SET12  |    SEA|       SEA89|             |             |      -0.1999|      ambiguous|
|SEA_SET12  |    SEA|       SEA90|             |             |       0.2877|      ambiguous|
|SWP_SET1   |    SWP|       SWP01|             |             |       0.6539|      ambiguous|
|SWP_SET1   |    SWP|       SWP02|             |             |       0.2333|      ambiguous|
|SWP_SET1   |    SWP|       SWP03|             |             |       0.5145|      ambiguous|
|SWP_SET1   |    SWP|       SWP04|             |             |       0.9011|      ambiguous|
|SWP_SET1   |    SWP|       SWP05|             |             |       0.7080|      ambiguous|
|SWP_SET1   |    SWP|       SWP06|       0.9597|       0.9202|       0.9644|      ambiguous|
|SWP_SET1   |    SWP|       SWP07|             |             |       0.6812|      ambiguous|
|SWP_SET1   |    SWP|       SWP08|             |             |             |      ambiguous|
|SWP_SET1   |    SWP|       SWP09|             |             |       0.7778|      ambiguous|
|SWP_SET2   |    SWP|       SWP10|       0.6619|       0.6875|       0.9078|      ambiguous|
|SWP_SET2   |    SWP|       SWP11|             |             |       0.5686|      ambiguous|
|SWP_SET2   |    SWP|       SWP12|       0.8882|       0.8370|       0.7923|      ambiguous|
|SWP_SET2   |    SWP|       SWP13|       0.4567|       0.7540|       0.3073|      ambiguous|
|SWP_SET2   |    SWP|       SWP14|             |             |       0.2852|      ambiguous|
|SWP_SET2   |    SWP|       SWP15|             |             |       0.8808|      ambiguous|
|SWP_SET2   |    SWP|       SWP16|       0.9503|       0.9064|       0.8841|      ambiguous|
|SWP_SET3   |    SWP|       SWP17|             |             |       0.8002|      ambiguous|
|SWP_SET3   |    SWP|       SWP18|             |             |       0.7583|      ambiguous|
|SWP_SET3   |    SWP|       SWP19|             |             |       0.2976|      ambiguous|
|SWP_SET3   |    SWP|       SWP20|       0.5537|       0.9771|       0.5913|      ambiguous|
|SWP_SET3   |    SWP|       SWP21|             |             |       0.5889|      ambiguous|
|SWP_SET3   |    SWP|       SWP22|             |             |       0.2520|      ambiguous|
|SWP_SET3   |    SWP|       SWP23|             |             |       0.7589|      ambiguous|
|SWP_SET3   |    SWP|       SWP24|             |             |             |      ambiguous|
|SWP_SET3   |    SWP|       SWP25|             |             |             |      ambiguous|
|SWP_SET3   |    SWP|       SWP26|             |             |       0.9470|      ambiguous|
|SWP_SET4   |    SWP|       SWP27|             |             |       0.8832|           cont|
|SWP_SET4   |    SWP|       SWP28|             |             |       0.6660|           cont|
|SWP_SET4   |    SWP|       SWP29|       0.3027|       0.9778|       0.3674|           cont|
|SWP_SET4   |    SWP|       SWP30|       0.6840|       0.9621|       0.7383|           cont|
|SWP_SET4   |    SWP|       SWP31|       0.6371|       0.8314|       0.8390|           cont|
|SWP_SET4   |    SWP|       SWP32|       0.6071|       0.9837|       0.6170|           cont|

Table: Cross-correlations between methodology pairs (Obs--Sim, Obs--Krg, Sim--Krg) for each individual SubMap transect. Empty cells indicate insufficient data to compute a correlation. NPA: North Pacific--Aleutians; SAM: South America; SEA: Southeast Asia; SWP: Southwest Pacific. {#tbl:method-ccf-summary}

\clearpage

|SubMap Set | Region| Model| Iteration| Max Distance| Vgram Cutoff| Vgram Lags| Observations| Max Pairs| Total Cost|
|:----------|------:|-----:|---------:|------------:|------------:|----------:|------------:|---------:|----------:|
|NPA_SET1   |    NPA|   Exp|       572|       131.90|          3.0|       20.0|         1170|       5.0|     0.9038|
|NPA_SET2   |    NPA|   Exp|       514|       373.80|          3.4|       26.0|          624|      47.9|     0.6761|
|NPA_SET3   |    NPA|   Exp|       507|       111.12|         11.3|       20.8|          339|       2.5|     0.9105|
|NPA_SET4   |    NPA|   Gau|       485|       105.48|          9.5|       30.8|          150|      31.5|     0.9155|
|NPA_SET5   |    NPA|   Exp|       464|       178.10|          4.3|       20.7|           94|      32.4|     0.8261|
|NPA_SET6   |    NPA|   Gau|       490|       104.59|          9.7|       20.3|           43|      50.0|     0.7040|
|NPA_SET7   |    NPA|   Sph|       492|       105.90|          6.4|       23.9|           79|      26.6|     0.8743|
|NPA_SET8   |    NPA|   Exp|       569|       439.92|         10.0|       20.6|         5972|      22.5|     0.7712|
|SAM_SET1   |    SAM|   Gau|       403|       309.95|          8.0|       39.5|          663|      24.8|     0.6510|
|SAM_SET2   |    SAM|   Sph|       516|       487.87|          9.4|       20.0|          855|      43.8|     0.6008|
|SAM_SET3   |    SAM|   Exp|       483|       246.00|          3.7|       22.8|          152|      44.0|     0.8659|
|SAM_SET4   |    SAM|   Gau|       377|       460.35|          4.3|       20.7|          631|      50.0|     0.8326|
|SAM_SET5   |    SAM|   Exp|       511|       425.95|          4.6|       23.0|          612|      21.9|     0.5713|
|SAM_SET6   |    SAM|   Exp|       474|       105.90|         11.5|       27.2|           94|      43.6|     0.7815|
|SAM_SET7   |    SAM|   Sph|       481|       156.39|         11.1|       21.4|          128|      49.9|     0.6717|
|SAM_SET8   |    SAM|   Gau|       493|       449.89|          3.2|       24.6|          169|      27.8|     0.7554|
|SAM_SET9   |    SAM|   Sph|       526|       119.49|          4.5|       47.2|          204|       9.6|     0.7751|
|SAM_SET10  |    SAM|   Sph|       505|       159.36|         10.5|       31.6|           73|       2.7|     0.8430|
|SEA_SET1   |    SEA|   Sph|       479|       387.76|         11.9|       21.0|          489|       8.2|     0.7133|
|SEA_SET2   |    SEA|   Exp|       525|       106.64|          4.0|       30.0|          519|      40.4|     0.6005|
|SEA_SET3   |    SEA|   Exp|       509|       116.71|          5.8|       20.1|          230|       4.6|     0.8953|
|SEA_SET4   |    SEA|   Sph|       487|       215.89|          7.9|       20.6|          126|      48.6|     0.6930|
|SEA_SET5   |    SEA|   Exp|       479|       103.13|          3.0|       49.0|          235|      37.0|     0.8231|
|SEA_SET6   |    SEA|   Exp|       492|       334.15|          4.8|       20.3|          505|      49.9|     0.7286|
|SEA_SET7   |    SEA|   Exp|       506|       122.54|          4.6|       29.7|          439|       2.6|     0.8977|
|SEA_SET8   |    SEA|   Sph|       513|       394.94|          3.0|       20.0|          620|       2.2|     0.9210|
|SEA_SET9   |    SEA|   Sph|       540|       499.86|          6.8|       20.0|         1571|       2.0|     0.8812|
|SEA_SET10  |    SEA|   Sph|       504|       360.88|          5.6|       20.2|          778|      12.7|     0.8587|
|SEA_SET11  |    SEA|   Gau|       488|       115.88|          6.3|       20.1|          169|      37.7|     0.9549|
|SEA_SET12  |    SEA|   Exp|       474|       109.71|          6.3|       27.1|           54|      32.6|     0.8718|
|SWP_SET1   |    SWP|   Exp|       489|       100.91|          3.3|       25.3|          188|      12.3|     0.6974|
|SWP_SET2   |    SWP|   Sph|       488|       103.66|          4.1|       20.9|          101|      48.9|     0.8564|
|SWP_SET3   |    SWP|   Gau|       500|       100.17|         11.0|       32.9|          169|      16.6|     0.9486|
|SWP_SET4   |    SWP|   Gau|       520|       373.69|          7.8|       21.0|          390|      18.9|     0.9900|

Table: Optimized Kriging parameters and variogram models for all 34 SubMap transect-set domains. Model abbreviations: Exp = exponential, Sph = spherical, Gau = Gaussian. Max Distance is the maximum neighborhood size $d_\mathrm{max}$ (km); Vgram Cutoff is the lag cutoff constant $\kappa$; Vgram Lags is the number of lag intervals $n$; Observations is the number of IHFC 2024 observations within the convex polygon domain; Max Pairs is the maximum neighborhood observations $m$; Total Cost is the composite optimization cost (dimensionless, range 0--1). {#tbl:nlopt-summary}

\clearpage

# References {.unnumbered #sec:references}

::: {#refs}
:::
