#!/bin/bash
#SBATCH --job-name=getsites
#SBATCH --nodelist=yournode
#SBATCH --ntasks=1
#SBATCH --array=1-23%10
#SBATCH --mail-user=youremail
#SBATCH --mail-type=END
python get_sites.py -g hg38.fa -c chr"$SLURM_ARRAY_TASK_ID" --merge --gzip -n GRCh38