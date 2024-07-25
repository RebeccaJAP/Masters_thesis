#!/bin/bash -l
#SBATCH --job-name=inf_karytyped
#SBATCH --account=project_2010414
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=70GB
#SBATCH --array=1-16
#SBATCH --output=out/output%a.txt
#SBATCH --error=err/errors%a.txt

# Load r-env
module load r-env

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2010414" >> ~/.Renviron

dataname=$(sed -n "$SLURM_ARRAY_TASK_ID"p ref_annotated_data.txt)

# Run the R script
srun apptainer_wrapper exec Rscript --no-save inf_karyotyped.R $dataname

