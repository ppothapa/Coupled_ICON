#/bin/bash

# Checks whether the input files contain the OpenMP conditional compilation
# sentinels (i.e. '!$'). The command line arguments can be paths either to
# regular files or to directories. In the latter case, the directories are
# searched recursively and all files with names that match known patterns are
# considered as input files for checking. If no arguments are provided, the
# script runs the standard Buildbot test, i.e. checks files in ICON source
# directories.

# The known patterns are (space-separated list of single-quoted patterns):
known_patterns="'*f90' '*.F90' '*.inc' '*.incf'"

# ICON source directories (space-separated list of single-quoted paths relative
# to the root repo directory):
icon_directories="'src' 'support' 'externals/jsbach' 'externals/emvorado' 'externals/art'"

# Number of parallel jobs:
job_num=8

set -eu
set -o pipefail

list_files()
{
  issue_warn=$1; shift
  for input in "$@"; do
    if test -d "${input}"; then
      eval "find "${input}" -type f -a \( $(
        eval "set dummy ${known_patterns}; shift"
        for pattern in "$@"; do
          echo -n "-name '${pattern}' -o "
        done
        echo "-false") \)"
    elif test -f "${input}"; then
      echo "${input}"
    elif test x"$issue_warn" = xyes; then
      echo "WARNING: input argument '${input}' is neither a directory not a file" >&2
    fi
  done
}

if test $# -eq 0; then
  icon_dir=$(unset CDPATH; cd "$(dirname "$0")/../../.."; pwd)
  eval "set dummy $(
    eval "set dummy ${icon_directories}; shift"
    for dir in "$@"; do
      echo -n "'${icon_dir}/${dir}' "
    done); shift"
fi

exitcode=0

list_files yes "$@" | xargs -P ${job_num} -I{} -- ${SHELL} -c 'grep --color="auto" -HnP "^\s*!\\\$\s" {} >&2; test $? -eq 1' || {
  {
    echo "ERROR: input files contain OpenMP conditional compilation sentinels (see above)"
    echo "       replace the sentinels with the macro '#ifdef _OPENMP/#endif' directives"
  } >&2
  exitcode=1
}

exit ${exitcode}
