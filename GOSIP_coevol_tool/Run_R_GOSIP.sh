#!/bin/bash
#SBATCH --job-name=Rtool
#SBATCH --mail-user=<useremail>
#SBATCH --mail-type=ALL
#SBATCH --account=account
#SBATCH -o ./logs/Report.Rtool_%a.out
#SBATCH --qos=qos
#SBATCH --output=out_Rtool_log_%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1  
#SBATCH --mem=128gb
#SBATCH --time=4-00:00:00
#SBATCH --partition=partition



ml R/4.2

Rscript coevol_simulator.R

