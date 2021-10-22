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
      unset itr
      echo '\nMax iterations for optimization search?'
      echo 'Note: computation cost is high!'
      echo 'Recommended < 50 iterations'
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

      unset ncores
      echo '\nPlease enter number of cores for parallel computing [0 for default]'
      echo 'Available cores on this machine:'
      nproc --all
      read 'ncores?Cores: '
      while [[ -z ${ncores} ]]; do
        read 'ncores?Cores: '
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
  ls -1q data/opt*.RData | xargs -n 1
  echo 'Run another round of kriging?'
  read 'p?yes/no: '
  while true; do
    if [[ $p == 'yes' ]]; then
      unset itr
      echo '\nMax iterations for optimization search?'
      echo 'Note: computation cost is high!'
      echo 'Recommended < 50 iterations'
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

      unset ncores
      echo '\nPlease enter number of cores for parallel computing [0 for default]'
      echo 'Available cores on this machine:'
      nproc --all
      read 'ncores?Cores: '
      while [[ -z ${ncores} ]]; do
        read 'ncores?Cores: '
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
      echo 'Please choose an opt file to run [type full path]:'
      ls -1q data/opt*.RData | xargs -n 1
      while [[ ! -f $fname ]]; do
        read 'fname?Filename: '
        if [[ ! -f $fname ]]; then
          echo 'File does not exist, please try again (or cntrl+c to escape)'
          ls -1q data/opt*.RData | xargs -n 1
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