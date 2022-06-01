#!/bin/bash
#=============================================================================
#SBATCH --job-name=build-icon
#SBATCH --nodes=1
#SBATCH --time=01:00:00

#set -x

SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)

if test ! -z "${SLURM_SUBMIT_DIR}"; then
  ICON_DIR="${SLURM_SUBMIT_DIR}"
else
  ICON_DIR=$(cd "${SCRIPT_DIR}/../.."; pwd)
fi
echo "SCRIPT_DIR:        |${SCRIPT_DIR}|"
echo "SLURM_SUBMIT_DIR:  |${SLURM_SUBMIT_DIR}|"
echo "ICON_DIR:          |${ICON_DIR}|"
echo "BB_NAME:           |${BB_NAME}|"
echo "BB_SYSTEM:         |${BB_SYSTEM}|"
echo "BB_experimentList: |${BB_EXPERIMENTLIST}|"
echo "BB_experiment:     |${BB_EXPERIMENT}|"


defaultExperimentList='icon-dev'

experimentList=${experimentList:-${BB_EXPERIMENTLIST:-${defaultExperimentList}}}
experimentLogfile="${ICON_DIR}/run/buildbotsExperiments.log"

echo "experimentList:    |${experimentList}|"
echo "experimentLogfile: |${experimentLogfile}|"


# plan: simplify the buildbot server config and take all decisions within the
# icon repo in a transparent way

# In case only an experiment is given {{{
if [ "x${BB_EXPERIMENT}" != "x" ]; then

  # first check if the given expriment file/template exists
  if [ ! -f "${ICON_DIR}/run/${BB_EXPERIMENT}}" ]; then
    echo "Cannot find experiment ${ICON_DIR}/run/${BB_EXPERIMENT}}"
    exit 1;
  fi

  # create a temporary list with the given experiments
  experimentList="buildbot_explist.tmp"
  # creata a new list with the given experiment
  cd ${ICON_DIR}/scripts/buildbot_scripts; ./create_all_builders "${experimentList}";
  cd ${ICON_DIR} ; ./scripts/buildbot_scripts/addexp --list "${experimentList}" "${BB_EXPERIMENT}"
# }}}
else # use the given/default experiment list  {{{

  # check if the given list is part of the repository. if not, then run the corresponding script to create the list
  cd ${ICON_DIR}/scripts/buildbot_scripts;
  if test ! -f experiment_lists/${experimentList}; then ./create_list_${experimentList} ; fi
fi #}}}

# list the expriments for this builder and list| this will be used by buildbot later on
cd ${ICON_DIR}; ./scripts/buildbot_scripts/lsexperiments --list ${experimentList} ${BB_NAME} > ${experimentLogfile}

cat ${experimentLogfile}

# build the model
cd ${ICON_DIR}; ./scripts/buildbot_scripts/build --list ${experimentList} ${BB_NAME}
