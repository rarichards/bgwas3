#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=8:mem=32gb

#set -o errexit
set -o pipefail
#set -o nounset
set -o xtrace
set -o errtrace

#PBS -V

PROJECT_DIR="/rds/general/user/grl2718/home/Projects/bgwas3"
module load anaconda3/personal
source activate ${PROJECT_DIR}/env


