#!/bin/bash
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH -p short
#SBATCH --mem=4000
#SBATCH -o logs/hostname_%j.out
#SBATCH -e logs/hostname_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@example.com
#SBATCH --array=1-100
export OMP_NUM_THREADS=4
uptime
btr-predict run_settings/settings_ROSMAP_LPOCV_BM9_46_Ordinal.json
uptime
