#!/bin/zsh
# Exit if any command fails
set -e

# Check for R dependencies
echo 'Checking for required R packages'
./packages.R

# Get opt files from data/opt*.RData
fnum=$(find data -name 'opt*' -print | wc -l)
if [[ ! $fnum -gt 0 ]]; then
  echo 'No opt files found'
  echo 'Run kriging?'
  read 'p?yes/no: '
  while true; do
    if [[ $p == 'yes' ]]; then
      echo 'Max iterations for optimization search?'
      echo 'Note: computation cost is high!'
      echo 'Recommended < 50 iterations'
      read 'itr?Number: '
      echo 'Please choose optimization algorithm:'
      echo 'for more info see: https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/'
      echo '1: NLOPT_GN_DIRECT_L   (Direct global no gradients)'
      echo '2: NLOPT_LN_SBPLX      (Local no gradients)'
      echo '3: NLOPT_LN_NELDERMEAD (Local no gradients)'
      echo '4: NLOPT_LN_BOBYQA     (Local no gradients)'
      echo '5: NLOPT_LN_COBYLA     (Local no gradients)'
      read 'alg?Choice: '
      echo 'Please enter weights for computing cost function'
      echo 'Interpolation and variogram weights must add to one'
      read 'iwt?Interpolation weight [0-1]: '
      read 'vwt?Variogram weight [0-1]: '
      echo 'Please enter number of cores for parallel computing [0 for default]'
      echo 'Available cores on this machine:'
      nproc --all
      read 'ncores?Cores: '
      echo 'Please enter number of folds for computing k-fold cross-validation'
      echo 'Note: computation cost is high!'
      echo 'Leave-one-out [default] is considerably slower than k-fold'
      echo 'Leave-one-out is more stable, however'
      echo 'Enter 0 for leave-one-out cross-validation [default]'
      read 'nfold?Number of folds [0=default]: '
      ./krige.R $itr $alg $iwt $vwt $ncores $nfold
      break
    elif [[ $p == 'no' ]]; then
      echo okay bye
      exit 0
    else
      read 'p?yes/no: '
    fi
  done
  echo 'Kriging successfull'
  echo 'Visualize results?'
  read 'p?yes/no: '
  fname=$(ls -1q data/opt*.RData | tail -n 1)
  while true; do
    if [[ $p == 'yes' ]]; then
      ./preprocess.R
      ./base_plots.R
      ./interpolation_plots.R $fname
      ./summary_plots.R $fname
      echo 'Finished!'
      exit 0
    elif [[ $p == 'no' ]]; then
      echo 'okay bye'
      exit 0
    else
      read 'p?yes/no: '
    fi
  done
else
  echo 'Found previous kriging results:'
  ls -1q data/opt*.RData | xargs -n 1 basename
  echo 'Run another round of kriging?'
  read 'p?yes/no: '
  while true; do
    if [[ $p == 'yes' ]]; then
      echo 'Max iterations for optimization search?'
      echo 'Note: computation cost is high!'
      echo 'Recommended < 50 iterations'
      read 'itr?Number: '
      echo 'Please choose optimization algorithm:'
      echo 'see https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/'
      echo 'for more information'
      echo '1: NLOPT_GN_DIRECT_L'
      echo '2: NLOPT_LN_SBPLX'
      echo '3: NLOPT_LN_NELDERMEAD'
      echo '4: NLOPT_LN_BOBYQA'
      echo '5: NLOPT_LN_COBYLA'
      read 'alg?Choice: '
      echo 'Please enter weights for computing cost function'
      echo 'Interpolation and variogram weights must add to one'
      read 'iwt?Interpolation weight [0-1]: '
      read 'vwt?Variogram weight [0-1]: '
      echo 'Please enter number of cores for parallel computing [0 for default]'
      echo 'Available cores on this machine:'
      nproc --all
      read 'ncores?Cores: '
      echo 'Please enter number of folds for computing k-fold cross-validation'
      echo 'Enter 0 for leave-one-out cross-validation [default]'
      read 'nfold?Number of folds [0=default]: '
      ./krige.R $itr $alg $iwt $vwt $ncores $nfold
      break
    elif [[ $p == 'no' ]]; then
      break
    else
      read 'p?yes/no: '
    fi
  done
  echo 'Visualize results?'
  read 'p?yes/no: '
  while true; do
    if [[ $p == 'yes' ]]; then
      echo 'Please choose an opt file to run:'
      ls -1q data/opt*.RData | xargs -n 1 basename
      while [[ ! -f data/$fname ]]; do
        read 'fname?Filename: '
        if [[ ! -f data/$fname ]]; then
          echo 'File does not exist, please try again (or cntrl+c to escape)'
          ls -1q data/opt*.RData | xargs -n 1 basename
        fi
      done
      ./preprocess.R
      ./base_plots.R
      ./interpolation_plots.R $fname
      ./summary_plots.R $fname
      echo 'Finished!'
      exit 0
    elif [[ $p == 'no' ]]; then
      echo 'okay bye'
      exit 0
    else
      read 'p?yes/no: '
    fi
  done
fi