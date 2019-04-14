#!/bin/sh

module purge
#. /etc/profile.d/modules.sh

module load matlab/2017a
export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXXXX)

export MATLABPATH=${MATLABPATH}:/${HOME}/${PROJECT}:${HOME}/MATLAB
#source ${HOME}/MATLAB/setpath.sh
export MATLABPATH="${MATLABPATH}:${HOME}/vbmc:${HOME}/vbmc/utils:${HOME}/bads"

#PROBLEMFOLDER="'${HOME}/neurobench-problems'"
PROBLEMFOLDER="[]"

#Check if running as an array job
if [[ ! -z "$PBS_ARRAYID" ]]; then
        IID=${PBS_ARRAYID}
fi
#Check if running as an array job
if [[ ! -z "$SGE_TASK_ID" ]]; then
        IID=${SGE_TASK_ID}
fi
if [[ ! -z "$SLURM_ARRAY_TASK_ID" ]]; then
        IID=${SLURM_ARRAY_TASK_ID}
fi

# Run the program

PARAMS=$(awk "NR==${IID} {print;exit}" ${INPUTFILE})

echo ${PARAMS} ${VERBOSE}

cat<<EOF | matlab -nodisplay
ibl_change_add2path
cd('${WORKDIR}');
example_mice={'$PARAMS'}
fit_all_unbiased
EOF
