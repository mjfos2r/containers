#!/bin/bash

# Run this to update the list of containers to build

find . -maxdepth 1 \( -path ".git" -o -path "." \) -prune -o -type d -print
