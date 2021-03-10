#!/usr/bin/env bash

theta=$1 # 1st argument is theta prefix
kTopo=$2 # 2nd argument is kTopo prefix
run=run_${theta}${kTopo%_} # run name
run_dir=$(pwd)/../runs/$run # full path to run directory

# Link shared input files to run directory
echo "Linking shared input"
shared_input=$(pwd)/../input/shared
for file in $(ls $shared_input); do
  ln -sf $shared_input/$file $run_dir/$file
done

# Link relevant generated input files to run directory. These are prefixed by
# either kTopoNN or kThetaNN.
echo "Linking generated input"
generated_input=$(pwd)/../input/generated
prefixes=(${kTopo} ${theta})
for prefix in ${prefixes[@]}; do
  for file in $(ls $generated_input | grep $prefix); do
    ln -sf $generated_input/$file $run_dir/${file#$prefix} # link and trim prefix
  done
done

# Set up a build directory
if [[ $3 == --build ]]; then
  # Create mod subdirectory if it doesn't exist
  mod_dir=$(pwd)/../build/mods/${theta%_}
  [ ! -d $mod_dir ] && mkdir -p $mod_dir

  # Link shared code to mod subdirectory
  echo "Linking shared code"
  shared_code_dir=$(pwd)/../code/shared
  for file in $(ls $shared_code_dir); do
    ln -sf $shared_code_dir/$file $mod_dir/$file
  done

  # Link generated code to run's code directory. I think we only need to wory
  # about changes in theta for this.
  echo "Linking generated code"
  generated_code_dir=$(pwd)/../code/generated
  prefix=${theta}
  for file in $(ls $generated_code_dir | grep $prefix); do
    ln -sf $generated_code_dir/$file $mod_dir/${file#$prefix}
  done
fi
