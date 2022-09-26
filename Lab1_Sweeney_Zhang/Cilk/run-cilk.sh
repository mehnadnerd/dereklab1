#!/bin/bash
# filename: run-cilk.sh

#SBATCH -J cmm                              # job name
#SBATCH -o %j                               # output and error file name (%j expands to jobID)
#SBATCH -n 2                                # total number of mpi tasks requested
#SBATCH -N 1                                # number of mpi nodes requested
#SBATCH -p normal                           # queue (partition) -- normal, development, etc.
#SBATCH -t 00:00:5                          # run time (hh:mm:ss) - 30 seconds
#SBATCH --mail-user=jiahanzhang@utexas.edu
#SBATCH --mail-type=begin                   # email me when the job starts
#SBATCH --mail-type=end                     # email me when the job finishes
CILK_NWORKERS=2 ./test_mm 1 0 32
