#!/bin/bash

#PBS -l nodes=1:ppn=20,walltime=6:45:01 -I -q gpu

module load openblas

cd  /home/rpanda/ising/
make -B

./montecarlo.x > output
