#!/bin/bash
#SBATCH --job-name=threading_test
#SBATCH --partition=owners
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=nticea@stanford.edu
#SBATCH --output=threading_test.txt
#SBATCH --error=threading_test.txt
#SBATCH --open-mode=append

# load Julia module
ml julia

# run the Julia application
julia ../../../code/ThreeBandHubbardPhonons.jl/examples/threading_test.jl
