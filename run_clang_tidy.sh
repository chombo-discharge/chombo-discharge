#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
#
# SPDX-License-Identifier: GPL-3.0-or-later

# Run clang-tidy in parallel.
# Usage:
#   ./run_clang_tidy.sh              -- run on all project files
#   ./run_clang_tidy.sh file1.cpp    -- run on specific files (pre-commit / CI)
set -euo pipefail

JOBS=$(nproc 2>/dev/null || echo 4)
BUILD_DIR=build
DIR_FILTER='^(Source|Exec|Physics|Geometries)/.*\.cpp$'

if [ $# -eq 0 ]; then
    exec run-clang-tidy -p "$BUILD_DIR" -j"$JOBS" \
        'Source/|Exec/|Physics/|Geometries/'
fi

cpp_files=$(printf '%s\n' "$@" | grep -E "$DIR_FILTER" || true)
if [ -z "$cpp_files" ]; then
    exit 0
fi
echo "$cpp_files" | xargs -P"$JOBS" -n1 clang-tidy -p "$BUILD_DIR"
