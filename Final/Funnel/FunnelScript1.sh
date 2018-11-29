#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=5:00:00
cd $PBS_O_WORKDIR
module add languages/R-3.4.1-ATLAS
R CMD BATCH --no-save --quiet --slave --no-restore FunnelNSMeth1.R