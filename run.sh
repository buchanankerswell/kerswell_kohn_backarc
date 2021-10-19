#!/bin/zsh
echo 'Compiling data ...'
Rscript preprocess.R
echo '\nPlotting global summary ...'
Rscript global_plots.R
echo '\nConstructing variograms ...'
Rscript construct_vgrms.R
echo '\nSummarizing interpolations ...'
Rscript summary_plots.R
echo '\nVisualizing results ...'
Rscript visualize_diff.R
Rscript visualize_diff_comp.R
echo '\nDone'
echo '\nRun ./knit.sh to render manuscript'