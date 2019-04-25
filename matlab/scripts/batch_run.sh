#!/bin/bash
PROJECT="ibl-changepoint"
SHORTNAME=IBL
BASEDIR="${HOME}/${PROJECT}/matlab"
SOURCEDIR="${BASEDIR}/"
JOBSCRIPT="${BASEDIR}/scripts/myjob.sh"

TRAINSET="endtrain"

#Job parameters
RUN=$TRAINSET
INPUTFILE="${BASEDIR}/scripts/joblist_${RUN}.txt"
MAXID=$(sed -n $= ${INPUTFILE})

RUNTIME=6:00:00
MAXRT=NaN
VERBOSE=0
MAXFUNMULT="[]"

NODES="1"
PPN="1"
MEM="4000MB"
RESOURCES="nodes=${NODES}:ppn=${PPN},mem=${MEM},walltime=${RUNTIME}"

JOBLIST="1-$MAXID"

#Convert from commas to spaces
JOBLIST=${JOBLIST//,/ }
echo JOBS $JOBLIST

WORKDIR="${SCRATCH}/${PROJECT}/${RUN}"
mkdir ${WORKDIR}
cd ${WORKDIR}

JOBNAME=${SHORTNAME}${RUN}

# running on Prince
sbatch --error=slurm-%A_%a.err --verbose --array=${JOBLIST} --mail-type=FAIL --mail-user=${USER}@nyu.edu --mem=${MEM} --time=${RUNTIME} --nodes=${NODES} --ntasks-per-node=${PPN} --export=PROJECT=${PROJECT},RUN=${RUN},TRAINSET=${TRAINSET},MAXID=$MAXID,WORKDIR=$WORKDIR,USER=$USER,MAXRT=$MAXRT,INPUTFILE=${INPUTFILE},VERBOSE=${VERBOSE},MAXFUNMULT=${MAXFUNMULT} --job-name=${JOBNAME} ${JOBSCRIPT}
