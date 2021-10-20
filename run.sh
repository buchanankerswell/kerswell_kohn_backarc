#!/bin/zsh

# Exit if command fails
set -e

if [[ ! `ls -1 data/opt*.RData | wc -l ` -gt 0 ]]; then
  echo 'No opt files found'
  echo 'Run kriging?'
  read 'p?yes/no: '
  a='return'
  while [[ ! $a == 'pass' ]]; do
    if [[ $p == 'yes' ]]; then
      a='pass'
      echo 'Max iterations for optimization search?'
      read 'itr?Number: '
      echo 'Please choose optimization algorithm:'
      echo 'for more info see: https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/'
      echo '1: NLOPT_GN_DIRECT_L'
      echo '2: NLOPT_LN_SBPLX'
      echo '3: NLOPT_LN_NELDERMEAD'
      echo '4: NLOPT_LN_BOBYQA'
      echo '5: NLOPT_LN_COBYLA'
      read 'alg?Choice: '
      ./krige.R $itr $alg
    elif [[ $p == 'no' ]]; then
      break
    else
      read 'p?yes/no: '
    fi
  done
  echo 'Kriging successfull'
  echo 'Visualize results?'
  read 'p?yes/no: '
  a='return'
  fname=$(ls -1q data/opt*.RData | tail -n 1)
  while [[ ! $a == 'pass' ]]; do
    if [[ $p == 'yes' ]]; then
      a='pass'
      ./preprocess.R
      ./base_plots.R
      ./interpolation_plots.R $fname
      ./summary_plots.R $fname
    elif [[ $p == 'no' ]]; then
      break
    else
      read 'p?yes/no: '
    fi
  done
else
  echo 'Found previous kriging results:'
  ls -1q data/opt*.RData | xargs -n 1 basename
  echo 'Run another round of kriging?'
  read 'p?yes/no: '
  a='return'
  while [[ ! $a == 'pass' ]]; do
    if [[ $p == 'yes' ]]; then
      a='pass'
      echo 'Max iterations for optimization search?'
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
      ./krige.R $itr $alg
    elif [[ $p == 'no' ]]; then
      break
    else
      read 'p?yes/no: '
    fi
  done
  echo 'Visualize results?'
  read 'p?yes/no: '
  a='return'
  while [[ ! $a == 'pass' ]]; do
    if [[ $p == 'yes' ]]; then
      a='pass'
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
    elif [[ $p == 'no' ]]; then
      break
    else
      read 'p?yes/no: '
    fi
  done
fi
echo 'Finished!'
#echo '\nRun ./knit.sh to render manuscript'