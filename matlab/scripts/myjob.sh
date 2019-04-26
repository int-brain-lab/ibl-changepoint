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

echo ${PARAMS} ${TRAINSET} ${VERBOSE}

cat<<EOF | matlab -nodisplay
ibl_changepoint_add2path
cd('${WORKDIR}');
runtype='${RUNTYPE}'

switch runtype
	case 'train'
		mice_list={'$PARAMS'}
		train_set='${TRAINSET}'
		fit_all_unbiased
	case 'fit'
		model_names = {'changepoint_doublenoise_loadnoiseendtrain_runlength_probs','changepoint_doublenoise_loadnoiseunbiased_runlength_probs'};
		Nopts = 10;
		vbmc_flag = true;
		refit_flag = true;
		modelfits = batch_model_fit(model_names,'$PARAMS',Nopts,vbmc_flag,refit_flag);
end

EOF
