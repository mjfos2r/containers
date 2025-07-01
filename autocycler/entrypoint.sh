#!/bin/bash

# entrypoint for autocycler docker container. Fixes the following:
#  - If this container is run as interactive, everything is groovy,
#    however...
#  - If this container is run non-interactively within a script,
#    the lack of TERM definition causes problems with some of the assemblers.
#    ...cough unicycler cough...
#
# This entrypoint solves that, it also retains the default behavior
# written by Ryan Wick in the original:
# https://github.com/rrwick/Autocycler/tree/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick

set -e

export TERM="${TERM:-dumb}"
export PATH="/opt/conda/envs/autocycler/bin:$PATH"
export PLASSEMBLER_DB="/usr/local/share/plassembler/db/"

# if no args specified, print help, list the installed tools, then exit.
if [[ $# -eq 0 ]]; then
    autocycler --help
    echo -e "\nInstalled Tools:"
    cat /tmp/installed_tools.txt || echo "None found!"
    echo -e "\nMissing Tools:"
    cat /tmp/missing_tools.txt || echo "None missing!"
    exit 0
fi

# otherwise execute the passed command

if [[ $# -eq 1 && "$1" == *" "* ]]; then
    exec bash -c "$1"
else
    exec "$@"
fi
