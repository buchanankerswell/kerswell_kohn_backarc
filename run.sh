#!/bin/zsh
echo 'Compiling data ...'
Rscript compile.R
echo '\nPlotting global summary ...'
Rscript global_plots.R
echo '\nInterpolating differences ...'
Rscript interp_diff.R
echo '\nSummarizing interpolations ...'
Rscript summary.R
echo '\nVisualizing results ...'
Rscript visualize_diff.R
Rscript visualize_diff_comp.R
echo '\nDone'
echo '\nRun ./knit.sh to render manuscript'