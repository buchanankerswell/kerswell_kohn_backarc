\clearpage

# Kriging System and Optimization {.unnumbered #sec:kriging-system-and-optimization}

## Ordinary Kriging

This study applies local isotropic ordinary Kriging methods under the following general assumptions:

- $\hat{\gamma}(h)$ is directionally invariant (isotropic)
- $\hat{\gamma}(h)$ is evaluated in two-dimensions and neglects elevation
- The first and second moments of $Z(u)$ are assumed to follow the conditions:

\begin{equation}
  \begin{aligned}
    &E[Z(u)] = \hat{Z}(u) = constant \\
    &E[(Z(u + h) - \hat{Z}(u))(Z(u) - \hat{Z}(u))] = C(h)
  \end{aligned}
  \label{eq:krige-assumptions}
\end{equation}

where $h$ is the lag distance, $C(h)$ is the covariance function, $E[Z(u)]$ is the expected value of the random variable $Z(u)$, and $\hat{Z}(u)$ is the arithmetic mean of $Z(u)$.

Equation \ref{eq:krige-assumptions} is known as "weak second-order stationarity". It assumes the underlying probability distribution of the observations $Z(u)$ does not change in space and the covariance $C(h)$ only depends on the distance $h$ between two observations. These assumptions are expected to be valid in cases where the underlying natural process is stochastic, spatially continuous, and has the property of additivity such that $\frac{1}{n}\sum_{i=1}^n Z(u_i)$ has the same meaning as $Z(u)$ [@bardossy1997].

The following are two illustrative cases where Equation \ref{eq:krige-assumptions} is likely valid:

> The thickness of a sedimentary unit with a homogeneous concentration of radioactive elements can be approximated by $q_s = q_b + \int A \,dz$, where $q_b$ is a constant heat flux entering the bottom of the layer and $A$ is the heat production within the layer with thickness $z$ [@furlong2013]. If one has two samples, $Z(u_1)$ = 31 mW/m$^2$ and $Z(u_2)$ = 30.5 mW/m$^2$, their corresponding thicknesses would be $Z'(u_1)$ = 1000 m and $Z'(u_2)$ = 500 m for $A$ = 0.001 mW/m$^3$ and $q_b$ = 30 mW/m$^2$. The variable, $Z(u)$, in this case is additive because the arithmetic mean of the samples is a good approximation of the average sedimentary layer thickness, $(Z(u_1) + Z(u_2)) /$ 2 = 750 m.

> The age of young oceanic lithosphere can be approximated by $q_s(t) = kT_b(\pi\kappa t)^{-1/2}$, where $q_s(t)$ is surface heat flow of a plate with age, $t$, $T_b$ is the temperature at the base of the plate, $k$ is thermal conductivity, and $\kappa = k/\rho C_p$ is thermal diffusivity [@stein1992]. Using reasonable values for $k$ = 3.138 W/mK, $\rho$ = 3330 kg/m$^3$, $C_p$ = 1171 J/kgK, $T_b$ = 1350 $^\circ$C, two samples, $Z(u_1)$ = 180 mW/m$^2$ and $Z(u_2)$ = 190 mW/m$^2$, would correspond to plates with ages of $Z'(u_1)$ = 10 Ma, and $Z'(u_2)$ = 9 Ma, respectively. Since $Z(u_1) + Z(u_2) /$ 2 = 185 mW/m$^2$ and $Z'(185~mW/m^2)$ = 9.5 Ma = $Z'(u_1) + Z'(u_2) /$ 2, the variable $Z(u)$ in this case is also additive.

Equation \ref{eq:krige-assumptions} is likely invalid in regions that transition among two or more tectonic regimes, however. For example, the expected (mean) heat flow $E[Z(u)]$ will change when moving from a spreading center to a subduction zone and thus $E[Z(u)] \neq constant$ over the region of interest. In other words, stationarity is violated and Kriging estimates may become spurious. Careful selection of Kriging parameters (outlined below; e.g. maximum point-pairs to use for local Kriging) can reduce or eliminate violations of stationarity assumptions embodied in \ref{eq:krige-assumptions}.

The second step is fitting a variogram model $\gamma(h)$ to the experimental variogram. This study fits six popular variogram models with sills (or theoretical sills) to the experimental variogram. The models are defined as [@pebesma2004]:

\begin{equation}
  \begin{aligned}
    Bes &\leftarrow \gamma(h) = 1 - \frac{h}{a}\ K_1\left(\frac{h}{a}\right) \quad \text{for } \  h \geq 0 \\
    Cir &\leftarrow \gamma(h) =
    \begin{cases}
      \frac{2}{\pi}\frac{h}{a}\ \sqrt{1-\left(\frac{h}{a}\right)^2} + \frac{2}{\pi}\ arcsin\left(\frac{h}{a}\right) \quad \text{for } \  0 \leq h \leq a \\
      nug + sill \quad \text{for } \  h > a
    \end{cases} \\
    Exp &\leftarrow \gamma(h) = 1 - exp\left(\frac{-h}{a}\right) \quad \text{for } \  h \geq 0 \\
    Gau &\leftarrow \gamma(h) = 1 - exp\left(\left[\frac{-h}{a}\right]^2\right) \quad \text{for } \  h \geq 0 \\
    Lin &\leftarrow \gamma(h) =
    \begin{cases}
      \frac{h}{a} \quad \text{for } \  0 \leq h \leq a \\
      nug + sill \quad \text{for } \  h > a
    \end{cases} \\
    Sph &\leftarrow \gamma(h) =
    \begin{cases}
      \frac{3}{2}\frac{h}{a} - \frac{1}{2}\left(\frac{h}{a}\right)^3 \quad \text{for } \  0 \leq h \leq a \\
      nug + sill \quad \text{for } \  h > a
    \end{cases} \\
  \end{aligned}
  \label{eq:variogram-models}
\end{equation}

where $h$ is the lag distance, $nug$ is the nugget, $sill$ is the sill, $a$ is the effective range, $K_1$ is a modified Bessel function. The models are Bessel, Circular, Exponential, Gaussian, Linear, and Spherical. For models without explicit sills (Bes, Exp, Gau), the effective range $a$ is the distance where the variogram reaches 95% of its maximum defined as 4$a$, 3$a$, and $\sqrt{3}a$ for Bes, Exp, and Gau, respectively [@graler2016; @pebesma2004]. The function `fit.variogram` in `gstat` is used to try all variogram models. The best model is selected by the minimum weighted least squares [@pebesma2004] error with weights proportional to the number of points in each lag divided by the squared lag distance $wt = N(h)_k/h_k^2$. Gaussian models produce spurious results in every case and are not included in the final analysis. Moreover, Circular models produce indistinguishable results from Spherical models, and so too were omitted from the final analysis.

Ordinary Kriging is used for interpolation, which estimates unknown observations $\hat{Z}(u)$ as a linear combination of all known observations [@bardossy1997]:

\begin{equation}
  \hat{Z}(u) = \sum_{i=1}^n \lambda_i Z(u_i)
  \label{eq:best-linear-unbiased-estimator}
\end{equation}

The conditions in Equation \ref{eq:krige-assumptions} set up a constrained minimization problem that can be solved with a system of linear equations. The expected value of $Z(u)$ is assumed to be the mean according to \ref{eq:krige-assumptions}, so the weights must be:

\begin{equation}
  \begin{aligned}
    E[\hat{Z}(u)] &= \sum_{i=1}^n \lambda_i E[Z(u_i)] \\
    \sum_{i=1}^n \lambda_i &= 1
  \end{aligned}
  \label{eq:unbiased}
\end{equation}

This constraint is known as the unbiased condition, which states that the sum of the weights must equal one. However, there is an infinite set of real numbers one could use for the weights, $\lambda_i$. The goal is to find the set of weights in Equation \ref{eq:best-linear-unbiased-estimator} that minimizes the estimation variance. This can be solved by minimizing the covariance function, $C(h)$ from Equation \ref{eq:krige-assumptions}:

\begin{equation}
  \begin{aligned}
    & \sigma^2(u) = Var[Z(u) - \hat{Z}(u)] = \\
    & E\left[(Z(u) - \sum_{i=1}^n \lambda_i Z(u_i))^2\right] = \\
    & E\left[Z(u)^2 + \sum_{j=1}^n \sum_{i=1}^n \lambda_j \lambda_i Z(u_j)Z(u_i) - 2 \sum_{i=1}^n \lambda_i Z(u_i)Z(u)\right] = \\
    & C(0) + \sum_{j=1}^n \sum_{i=1}^n \lambda_j \lambda_i C(u_i - u_j) - 2 \sum_{i=1}^n \lambda_i C(u_i - u)
  \end{aligned}
  \label{eq:minimize-variance}
\end{equation}

Minimizing Equation \ref{eq:minimize-variance} with respect to the unbiased condition (Equation \ref{eq:unbiased}), yields the best linear unbiased estimator [BLUE, @bardossy1997] for Equation \ref{eq:best-linear-unbiased-estimator} and together comprise the Kriging system of equations. The functions `krige` and `krige.cv` in `gstat` are used for surface heat flow interpolation and error estimation by k-fold cross-validation [@pebesma2004].

## Optimization with `nloptr`

Achieving accurate Kriging results depends on one's choice of many Kriging parameters, $\Theta$. In this study, we investigate a set of parameters:

\begin{equation}
  \Theta = \{model,\ n_{lag},\ cut,\ n_{max},\ shift\}
  \label{eq:params}
\end{equation}

where $model$ is one of the variogram models defined in Equation \ref{eq:variogram-models}, $n_{lag}$ is the number of lags, $cut$ is a lag cutoff proportionality constant, $n_{max}$ is the maximum point-pairs for local Kriging, and $shift$ is a horizontal lag shift constant. The lag cutoff constant $cut$ controls the maximum separation distance between pairs of points used to calculate the experimental variogram (i.e. the x-axis range or "width" of the experimental variogram). The horizontal lag shift constant $shift$ removes the first few lags from being evaluated by effectively shifting all lags to the left proportionally by $shift$. This is necessary to avoid negative ranges when fitting experimental variograms with anomalously high variances at small lag distances.

The goal is to find $\Theta$ such that the Kriging function $f(x_i;\ \Theta)$ gives the minimum error defined by a cost function $C(\Theta)$, which represents the overall goodness of fit of the interpolation. This study defines a cost function that simultaneously considers errors between the experimental variogram $\hat{\gamma}(h)$ and modelled variogram $\gamma(h)$, and between surface heat flow observations $Z(u_i)$ and Kriging estimates $\hat{Z}(u)$ [after @li2018]:

\begin{equation}
  \begin{aligned}
    C(\Theta) &= w_{vgrm}\ C_{vgrm}(\Theta) + w_{interp}\ C_{interp}(\Theta) \\
    &w_{vgrm} + w_{interp} = 1
  \end{aligned}
  \label{eq:cost}
\end{equation}

where $C_{vgrm}(\Theta)$ is the normalized RMSE evaluated during variogram fitting and $C_{interp}(\Theta)$ is the normalized RMSE evaluated during Kriging. Weighted ordinary least squares is used to evaluate $C_{vgrm}(\Theta)$, whereas k-fold cross-validation is used to evaluate $C_{interp}(\Theta)$. K-fold splits the dataset $|Z(u_i)|$ into $k$ equal intervals, removes observations from an interval, and then estimates the removed observations by fitting a variogram model to data in the remaining $k-1$ intervals. This process is repeated over all $k$ intervals so that the whole dataset has been cross-validated. The final expression to minimize becomes:

\begin{equation}
  \begin{aligned}
    &C(\Theta) = \\
    &\frac{w_{vgrm}}{\sigma_{vgrm}}\ \left(\frac{1}{N(h)}\ \sum_{k=1}^{N}\ w(h_k)\ [\hat{\gamma}(h_k)-\gamma(h_k;\ \Theta)]^2\right)^{1/4} + \\
    &\frac{w_{interp}}{\sigma_{interp}}\ \left(\frac{1}{M}\ \sum_{i=1}^{M}\ [Z(u_i)-\hat{Z}(u_i;\ \Theta)]^2\right)^{1/2}
  \end{aligned}
  \label{eq:cost-function}
\end{equation}

where $N(h)$ is the number of point-pairs used to evaluate the experimental variogram and $w(h_k) = N(h)_k/h_k^2$ are weights defining the importance of the $kth$ lag on the variogram model fit.  $Z(u_i)$ and $\hat{Z}(u_i;\ \Theta)$ are the observed and estimated values, respectively, and m is the number of measurements in $|Z(u_i)|$. The RMSEs are normalized by dividing by $\sigma_{vgrm}$ and $\sigma_{interp}$, which represent the standard deviation of the experimental variogram $\hat{\gamma}(h)$ and surface heat flow observations $Z(u_i)$, respectively. The weights $w_{vgrm}$ and $w_{interp}$ were varied between 0 and 1 to test the effects on $C(\Theta)$. Preferred weights of $w_{vgrm}$ = $w_{interp}$ = 0.5 are selected to balance the effects of $C_{vgrm}(\Theta)$ and $C_{interp}(\Theta)$ on the cost function.

![Summary of optimized Kriging parameters. Cost does not correlate strongly with most Kriging parameters (solid black line with ivory 95% confidence intervals), indicating the optimization procedure is successfully generalizable across subduction zone segments. The exception is a correlation between cost and the logarithm of the experimental variogram sill. Note that parameter values adjust from an initial value (solid white line) during the optimization procedure.](../figs/vgrmSummary.png){#fig:variogram-summary width=100%}

Minimization of $C(\Theta)$ is achieved by non-linear constrained optimization using algorithms defined in the R package `nloptr` [@ypma2014]. Global search methods had limited success compared to local search methods. See [the official documentation](https://nlopt.readthedocs.io/en/latest/NLopt_Introduction/) for more information on `nloptr` and available optimization algorithms. The run used to produce the visualizations in this study apply the `NLOPT_LN_COBYLA` method [constrained optimization by linear approximation, @powell1994] with 50 max iterations, leave-one-out cross-validation (k-fold $=$ the number of observations) in the evaluated segment, and cost function weights of $w_{vgrm}$ = $w_{interp}$ = 0.5 (Figure \ref{fig:nlopt-loss-curve}). All data, code, and instructions to reproduce results in this study can be found at [https://github.com/buchanankerswell/kerswell_kohn_backarc](https://github.com/buchanankerswell/kerswell_kohn_backarc).

![Cost function minimization for Kriging interpolations. Most variogram models (panels) converge on a local optimum for most Kriging domains (lines) after 15-20 iterations. Each line represents one of thirteen subduction zone segments. See text for bound constraints and other options passed to the optimization procedure.](../figs/nlopt-loss-curve.png){#fig:nlopt-loss-curve width=100%}

# Variogram Models

**!!!fitted-variogram-models.png**
**!!!variogram-summary-table.md**

# ThermoGlobe Summary

![Distribution of ThermoGlobe observations from @lucazeau2019 cropped within 1000 km-radius buffers around 13 active subduction zone segments. Heat flow distributions are centered between **!!!number** and **!!!number** mW/m$^2$, generally right-skewed, and irregularly distributed. Skewness reflects near-surface perturbations from geothermal systems and tectonic regions with high thermal activity while irregularity reflects complex heat exchange acting across multiple spatial scales from 10$^-1$ to 10$^3$ km.](../figs/hfSummary.png){#fig:heatflow-observations-summary width=100%}

**!!!heatflow-observations-summary-table.md**

# Comparing Similarity and Kriging Interpolations

![Differences between Similarity and Kriging interpolations by segment, computed as Similarity-Kriging. Differences are centered near zero with medians ranging from **!!!number** to **!!!number** mW/m$^2$, but broadly distributed with IQRs from **!!!number** to **!!!number** mW/m$^2$ and some long tails extending from **!!!number** to **!!!number** mW/m$^2$. Positive medians and right skew indicate a general tendency towards higher surface heat flow predictions by Similarity compared to Kriging. The broadest distributions (Andes and Central America) reflect less subtle differences between methods. Distributions are colored by quartiles (25%, 50%, 75%). Similarity interpolation from @lucazeau2019.](../figs/interpDiffSummary.png){#fig:interpolation-differences-summary width=100%}

![Summary of differences between Similarity and Kriging uncertainties computed as Similarity-Kriging. Differences are centered at slightly negative values with median differences ranging from **!!!number** to **!!!number** mW/m$^2$, and relatively narrowly distributed with IQRs from **!!!number** to **!!!number** mW/m$^2$ and some long tails extending from **!!!number** to **!!!number** mW/m$^2$. Negative medians indicate greater uncertainties by Kriging compared to Similarity. Distributions are colored by quantiles (25%, 50%, 75%). Similarity data from @lucazeau2019. Refer to Figure \ref{fig:interpolation-differences-summary} for estimate differences.](../figs/interpSigmaDiffSummary.png){#fig:sigmaDiffSummaryPlot width=100%}

**!!!diff-heatflow-summary-table.md**
**!!!diff-sigma-summary-table.md**
**!!!diff-segment-examples.png**

# Upper-plate Surface Heat Flow

**!!!heatflow-profiles-segment-examples.png**
**!!!sector-summary-table.md**

# References {.unnumbered #sec:references}

::: {#refs}
:::
