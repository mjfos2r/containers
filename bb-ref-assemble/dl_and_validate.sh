#!/bin/bash

set -euo pipefail

# Colors for output
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[0;33m'; NC='\033[0m'

verify_checksum() {
    echo -e "${YELLOW}Validating md5 checksum${NC}"
    local expected actual
    expected=$(cut -d' ' -f1 "$1")
    file=$2
    filename=$(basename "$file")

    if command -v md5sum &>/dev/null; then actual=$(md5sum "$file" | awk '{print $1}')
    else echo -e "[ ${RED}FAIL${NC} ]::[ HASH ]::[ ${RED}No md5 utility found!${NC} ]"; return 0;
    fi

    if [[ $expected == "$actual" ]]; then echo -e "[ ${GREEN}PASS${NC} ]::[ HASH ]::[ ${YELLOW}${filename}${NC} ]"
    else { echo -e "[ ${RED}FAIL${NC} ]::[ HASH ]::[ ${YELLOW}${filename}${NC} ]"; exit 1; }
    fi
}

BINARIES=(
    ""
)

