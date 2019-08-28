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

echo ${PARAMS} ${TRAINSET} ${SECONDPARAM}

cat<<EOF | matlab -nodisplay
ibl_changepoint_add2path
cd('${WORKDIR}');
runtype='${RUNTYPE}'

switch runtype
	case 'psycho'
		mice_list={'$PARAMS'}
		fit_psycho_sessions
	case 'train'
		mice_list={'$PARAMS'}
		train_set='${TRAINSET}'
		fit_all_unbiased
	case 'fit'
		% model_names = {'changepoint_contrastnoise_runlength_probs'};
		% model_names = {'changepoint_contrastnoise'};
		% model_names = {'omniscient_contrastnoise_fixedprior'};
		% model_names = {'exponential_contrastnoise'};
		model_names = {'$SECONDPARAM'};
		mouse_name = ['$PARAMS'];
		% mouse_name = ['$PARAMS' '_half2'];
		Nopts = [10,5];
        	hmm_flag = false;
		vbmc_flag = false;
		refit_flags = [false,false,false];
		modelfits = batch_model_fit(model_names,mouse_name,Nopts,hmm_flag,vbmc_flag,refit_flags);
	case 'learn'
		mouse_name = ['$PARAMS'];
		data_mod='${SECONDPARAM}';
		if data_mod(1) == '['; data_mod = []; end % Passed empty string
		run_changepoint_learner
end
EOF
