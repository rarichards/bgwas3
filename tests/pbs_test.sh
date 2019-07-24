#PBS -l walltime=00:01:00
#PBS -l select=1:ncpus=1:mem=1gb

#set -o errexit
set -o pipefail
#set -o nounset
set -o xtrace
set -o errtrace

#PBS -V

echo "pbs working"
