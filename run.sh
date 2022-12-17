#!/bin/zsh

# Clock time
SECONDS=0

# Exit if any command fails
set -e
# Check for R dependencies
R/packages.R
# Check for files in data directory
fnum=$(find data -name '*' -print | wc -l)
if [[ $fnum -lt 33 ]]; then
  # Download data from osf
  R/download-data.R
fi
# Get opt files from data/opt*.RData
fnum=$(find data -name 'hf.RData' -print | wc -l)
if [[ ! $fnum -gt 0 ]]; then
  # Preprocess data
  R/preprocess.R
  cp data/hf.RData draft/assets/r/
fi
# Get opt files from data/opt*.RData
fnum=$(find data -name 'opt*' -print | wc -l)
if [[ ! $fnum -gt 0 ]]; then
  echo '\nNo opt file found'
  echo 'Run kriging?'
  read 'p?yes/no: '
  while true; do
    if [[ $p == 'yes' ]]; then
      unset itr
      echo '\nMax iterations for optimization search?'
      echo 'Note: computation cost is high!'
      echo 'Recommended <= 30 iterations'
      read 'itr?Number: '
      while [[ -z ${itr} ]]; do
        read 'itr?Number: '
      done
      unset alg
      echo '\nPlease choose optimization algorithm:'
      echo 'for more info on algorithms see:'
      echo 'https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/'
      echo '\n1: NLOPT_GN_DIRECT_L   (Direct global no gradients)'
      echo '2: NLOPT_LN_SBPLX      (Local no gradients)'
      echo '3: NLOPT_LN_NELDERMEAD (Local no gradients)'
      echo '4: NLOPT_LN_BOBYQA     (Local no gradients)'
      echo '5: NLOPT_LN_COBYLA     (Local no gradients)'
      read 'alg?Choice: '
      while [[ -z ${alg} ]]; do
        read 'alg?Choice: '
      done
      unset iwt
      unset vwt
      echo '\nPlease enter weights for computing cost function'
      echo 'Note: interpolation and variogram weights must add to one'
      read 'iwt?Interpolation weight [0-1]: '
      while [[ -z ${iwt} ]]; do
        read 'iwt?Interpolation weight [0-1]: '
      done
      read 'vwt?Variogram weight [0-1]: '
      while [[ -z ${vwt} ]]; do
        read 'vwt?Variogram weight [0-1]: '
      done
      unset nfold
      echo '\nPlease enter number of folds for computing k-fold cross-validation'
      echo 'Note: computation cost is high!'
      echo 'Leave-one-out [default] is considerably slower than k-fold, but more stable'
      echo '\nEnter 0 for leave-one-out cross-validation [default]'
      echo 'or enter a value between 0 and 1 to set k-folds as a'
      echo 'proportion of the number of data points in the kriging domain'
      read 'nfold?Number of folds [0=default]: '
      while [[ -z ${nfold} ]]; do
        read 'nfold?Number of folds [0=default]: '
      done
      unset ncores
      echo '\nPlease enter number of cores for parallel computing [0 for default]'
      echo 'Available cores on this machine:'
      nproc --all
      read 'ncores?Cores: '
      while [[ -z ${ncores} ]]; do
        read 'ncores?Cores: '
      done
      R/krige.R $itr $alg $iwt $vwt $nfold $ncores
      cp data/opt.RData draft/assets/r/
      echo 'Kriging successfull!'
      break
    elif [[ $p == 'no' ]]; then
      echo okay bye
      # Print clock time
      t=$SECONDS
      printf '\nTime taken: %d days, %d minutes, %d seconds\n' \
        "$(( t/86400 ))" "$(( t/60 - 1440*(t/86400) ))" "$(( t ))"
      exit 0
    else
      read 'p?yes/no: '
    fi
  done
  echo '\nvisualize results?'
  read 'p?yes/no: '
  while true; do
    if [[ $p == 'yes' ]]; then
      R/base-plots.R
      R/interpolation-plots.R
      R/summary-plots.R
      R/goutorbe-analysis.R
      cp data/sectors.RData draft/assets/r/
      echo 'finished!'
      # print clock time
      t=$SECONDS
      printf '\ntime taken: %d days, %d minutes, %d seconds\n' \
        "$(( t/86400 ))" "$(( t/60 - 1440*(t/86400) ))" "$(( t ))"
      exit 0
    elif [[ $p == 'no' ]]; then
      echo 'okay bye'
      # print clock time
      t=$SECONDS
      printf '\ntime taken: %d days, %d minutes, %d seconds\n' \
        "$(( t/86400 ))" "$(( t/60 - 1440*(t/86400) ))" "$(( t ))"
      exit 0
    else
      read 'p?yes/no: '
    fi
  done
elif [[ $fnum -gt 0 ]]; then
  echo '\nopt file found at data/opt.RData'
  echo 'visualize results?'
  read 'p?yes/no: '
  while true; do
    if [[ $p == 'yes' ]]; then
      R/base-plots.R
      R/interpolation-plots.R
      R/summary-plots.R
      R/goutorbe-analysis.R
      cp data/sectors.RData draft/assets/r/
      echo 'finished!'
      # print clock time
      t=$SECONDS
      printf '\ntime taken: %d days, %d minutes, %d seconds\n' \
        "$(( t/86400 ))" "$(( t/60 - 1440*(t/86400) ))" "$(( t ))"
      exit 0
    elif [[ $p == 'no' ]]; then
      echo 'okay bye'
      # print clock time
      t=$SECONDS
      printf '\ntime taken: %d days, %d minutes, %d seconds\n' \
        "$(( t/86400 ))" "$(( t/60 - 1440*(t/86400) ))" "$(( t ))"
      exit 0
    else
      read 'p?yes/no: '
    fi
  done
fi