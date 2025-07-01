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

# force a dumb terminal if we're not running as interactive.
if [[ -z "$TERM" ]]; then
    export TERM=dumb
fi

tools=(
    "canu" "flye" "lja" "metamdbg"
    "miniasm" "minimap2" "minipolish"
    "necat" "nextdenovo" "nextpolish"
    "racon" "raven-assembler" "wtdbg"
    "plassembler" "unicycler"
)

# ensure the path contains the conda binaries.
export PATH="/opt/conda/envs/autocycler/bin:$PATH"
# and that plassembler_db is set.
export PLASSEMBLER_DB="/usr/local/share/plassembler/"

if [[ $# -eq 0 ]]; then
    autocycler --help
    echo -e "\nInstalled Tools:"
    cat /tmp/installed_tools.txt || echo "None found!"
    echo -e "\nMissing Tools:"
    cat /tmp/missing_tools.txt || echo "None missing!"
    echo -e"\nTool Versions:"
    for tool in "${tools[@]}"; do
        if command -v "$tool" &>/dev/null; then
            version=$("$tool" --version 2>&1 | head -n1)
            echo "[ PASS ]::[ $tool ]::[ $version ]"
        else
            echo "[ FAIL ]::[ $tool ]::[ not found! ]"
        fi
    done
else
    # otherwise execute the passed command
    exec "$@"
fi