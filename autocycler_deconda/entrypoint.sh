#!/usr/bin/env bash
# Keep it minimal so it can't interfere with Cromwell's wrapper

# no more set -e

export TERM="${TERM:-dumb}"
export PATH="/opt/conda/envs/autocycler/bin:$PATH"
export PLASSEMBLER_DB="/usr/local/share/plassembler/db/"

# If no args, show help and exit (still fine without set -e)
if [[ $# -eq 0 ]]; then
  autocycler --help || true
  echo -e "\nInstalled Tools:";  cat /tmp/installed_tools.txt 2>/dev/null || echo "None found!"
  echo -e "\nMissing Tools:";    cat /tmp/missing_tools.txt 2>/dev/null || echo "None missing!"
  exit 0
fi

# Always exec the given command/args exactly
exec "$@"
