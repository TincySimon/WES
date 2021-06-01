# Join stderr and stdout log files into stdout log file
#$ -j y
# Keep current environment variables
#$ -V
# Use current working directory as working directory of the job.  This
# requires you to cd/popd into the directory where the snake_job.sh lies.
#$ -cwd
# allocating end time
#$ -l h_rt=48:00:00
#$ -pe smp 6
#$ -l h_vmem=30g
#$ -m eas
#$ -o sge_log

# Enforce existence of TMPDIR -----------------------------------------------

export TMPDIR=/fast/users/${USER}/scratch/tmp
mkdir -p ${TMPDIR}
export _JAVA_OPTIONS=-Djava.io.tmpdir=${TMPDIR}
set -x

# Configuration variables Create one log directory per Snakemake run---------

LOGDIR=sge_log/${JOB_ID}
mkdir -p ${LOGDIR}

# activating miniconda3 -----------------------------------------------

source $HOME/scratch/miniconda3/bin/activate py27

# script  --------------------------------------------------------

/fast/projects/ngs_sers/2017_WES_Analysis/PNET/WES_workflow/Scripts/PublicationScripts/ScarpaAnalysis/ScarpaMutation_Annovar.sh
