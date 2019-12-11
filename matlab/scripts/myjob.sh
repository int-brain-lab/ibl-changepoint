#!/bin/sh

module purge
#. /etc/profile.d/modules.sh

module load matlab/2017a
export MATLAB_PREFDIR=$(mktemp -d $SLURM_JOBTMP/matlab-XXXXXX)

export MATLABPATH=${MATLABPATH}:/${HOME}/${PROJECT}:${HOME}/MATLAB
#source ${HOME}/MATLAB/setpath.sh
export MATLABPATH="${MATLABPATH}:${HOME}/vbmc:${HOME}/vbmc/utils:${HOME}/bads:${HOME}/hmmfit:${HOME}/vbmc-dev/acq:${HOME}/vbmc-dev/misc"

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
		% model_names = {'omniscient_contrastnoise_fixedfreeprior'};
		% model_names = {'exponential_contrastnoise'};
		model_names = {'$SECONDPARAM'};
		mouse_name = ['$PARAMS'];
		% mouse_name = ['$PARAMS' '_half2'];
		Nopts = [10,4];
        	hmm_flag = false;  vbmc_flag = true;  mcmc_flag = false;
        	methods_flags = [hmm_flag,vbmc_flag,mcmc_flag];
		refit_flags = [false,false,false,false];
		modelfits = batch_model_fit(model_names,mouse_name,Nopts,methods_flags,refit_flags);
    case {'prior1','prior2'}
		model_names = {'exponential_contrastnoise_hyperprobs'};
		mouse_name = ['$PARAMS'];
		Nopts = [10,4];
        	hmm_flag = false;  vbmc_flag = false;  mcmc_flag = true;
        	methods_flags = [hmm_flag,vbmc_flag,mcmc_flag];
		refit_flags = [false,false,false,false];
        	if strcmpi(runtype,'prior1'); suffix = 'endtrain'; else; suffix = 'unbiased'; end
        	modelfits = batch_model_fit(model_names,[mouse_name '_' suffix],Nopts,methods_flags,refit_flags);
	case 'learn'
		mouse_name = ['$PARAMS'];
		data_mod='${SECONDPARAM}';
		if data_mod(1) == '['; data_mod = []; end % Passed empty string
		run_changepoint_learner
end
EOF
