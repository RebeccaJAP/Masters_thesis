#!/bin/bash -l
#SBATCH --job-name=sce_karyotyped
#SBATCH --account=project_2010414
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=70GB
#SBATCH --array=1-35
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

# Change ref_annotated_data_clust.txt to ref_annotated_data_d20.txt to use the original cell types data
# and change the array size accordingly.
dataname=$(sed -n "$SLURM_ARRAY_TASK_ID"p ref_annotated_data_clust.txt)

# Run the R script
srun apptainer_wrapper exec Rscript --no-save sce_karyotyped.R $dataname

