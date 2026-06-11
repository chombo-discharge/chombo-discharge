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
DIR_FILTER='^(Source|Exec|Physics|Geometries)/.*\.cpp$'

# Rewrite compile_commands.json so every -I path that goes into Submodules
# becomes -isystem instead.  The compiler and clang-analyzer both suppress
# diagnostics from system headers, preventing third-party warnings from
# polluting the output.  Regular clang-tidy checks are already handled by
# HeaderFilterRegex in .clang-tidy; -isystem is needed specifically for
# clang-analyzer findings whose reported location lands inside a Submodule
# header (those bypass the header filter).
make_clean_db() {
    python3 - "$1" <<'PYEOF'
import json, re, sys
out = sys.argv[1]
with open('compile_commands.json') as f:
    cmds = json.load(f)
for c in cmds:
    if 'command' in c:
        c['command'] = re.sub(r'-I(/[^ ]*Submodules[^ ]*)', r'-isystem \1', c['command'])
    elif 'arguments' in c:
        new_args = []
        for a in c['arguments']:
            m = re.match(r'^-I(/.*Submodules.*)$', a)
            if m:
                new_args.extend(['-isystem', m.group(1)])
            else:
                new_args.append(a)
        c['arguments'] = new_args
with open(out + '/compile_commands.json', 'w') as f:
    json.dump(cmds, f)
PYEOF
}

# Use clang-tidy-cache as a transparent wrapper when available.
if command -v clang-tidy-cache &>/dev/null; then
    CLANG_TIDY_BIN="clang-tidy-cache"
else
    CLANG_TIDY_BIN="clang-tidy"
fi

if [ $# -eq 0 ]; then
    WORK_DIR=$(mktemp -d)
    trap "rm -rf $WORK_DIR" EXIT
    make_clean_db "$WORK_DIR"
    run-clang-tidy -p "$WORK_DIR" -j"$JOBS" \
        -clang-tidy-binary "$CLANG_TIDY_BIN" \
        'Source/|Exec/|Physics/|Geometries/'
    exit $?
fi

cpp_files=$(printf '%s\n' "$@" | grep -E "$DIR_FILTER" || true)
if [ -z "$cpp_files" ]; then
    exit 0
fi

WORK_DIR=$(mktemp -d)
trap "rm -rf $WORK_DIR" EXIT
make_clean_db "$WORK_DIR"
echo "$cpp_files" | xargs -P"$JOBS" -n1 "$CLANG_TIDY_BIN" -p "$WORK_DIR"
