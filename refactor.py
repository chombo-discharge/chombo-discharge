#!/usr/bin/env python3
import os
import re
import sys

# 1) Replace MPI_Init block with ChomboDischarge::initialize(argc, argv)
INIT_BLOCK_RE = re.compile(
    r'#ifdef\s+CH_MPI\s*\n'
    r'\s*MPI_Init\s*\(\s*&argc\s*,\s*&argv\s*\)\s*;\s*\n'
    r'\s*#endif',
    re.MULTILINE
)

# 2) Replace any CH_MPI block that contains MPI_Finalize() with ChomboDischarge::finalize();
FINALIZE_BLOCK_RE = re.compile(
    r'#ifdef\s+CH_MPI\b.*?#endif',
    re.DOTALL | re.MULTILINE
)

# 3a) Remove const std::string input_file = argv[1];
INPUT_FILE_RE = re.compile(
    r'^\s*const\s+std::string\s+input_file\s*=\s*argv\[\s*1\s*\]\s*;\s*$',
    re.MULTILINE
)

# 3b) Catch ANY ParmParse variable declaration (possibly multi-line)
#     Example: ParmParse pp(argc - 2, argv + 2, NULL, input_file.c_str());
PARMPARSE_DECL_RE = re.compile(
    r'\bParmParse\s+([A-Za-z_]\w*)\s*\((.*?)\)\s*;',
    re.DOTALL
)

# 4) Replace setupAndRun(<anything>) with setupAndRun()
SETUP_AND_RUN_RE = re.compile(
    r'setupAndRun\s*\([^)]*\)'
)

# 5) RefCountedPtr<Class> var = ... -> auto var = ...
REFCOUNTEDPTR_LHS_RE = re.compile(
    r'\bRefCountedPtr\s*<[^>]+>\s+([A-Za-z_]\w*)\s*='
)


def remove_parmparse_and_uses(code: str) -> str:
    """
    - Find all ParmParse variable declarations with 4 arguments
    - Remove the declarations
    - Remove all lines that use those variables
    """
    to_kill_vars = set()

    def decl_replacer(m: re.Match) -> str:
        var_name = m.group(1)
        args = m.group(2)
        # crude but effective: 4 args => at least 3 commas
        if args.count(',') >= 3:
            to_kill_vars.add(var_name)
            return ''   # drop the whole 'ParmParse var(...);' statement
        return m.group(0)

    # Remove matching ParmParse declarations and collect the variable names
    code = PARMPARSE_DECL_RE.sub(decl_replacer, code)

    if not to_kill_vars:
        return code

    # Remove any lines that reference those variables (their "uses")
    lines = code.splitlines()
    new_lines = []
    for line in lines:
        if any(re.search(r'\b' + re.escape(v) + r'\b', line) for v in to_kill_vars):
            # Drop this line; it's using a ParmParse we removed
            continue
        new_lines.append(line)

    return "\n".join(new_lines)


def transform_content(code: str) -> str:
    # --- 1) MPI_Init block -> ChomboDischarge::initialize(argc, argv);
    code = INIT_BLOCK_RE.sub('ChomboDischarge::initialize(argc, argv);', code)

    # --- 2) CH_MPI block with MPI_Finalize() -> ChomboDischarge::finalize();
    def replace_finalize_block(match: re.Match) -> str:
        block = match.group(0)
        if 'MPI_Finalize' in block:
            return 'ChomboDischarge::finalize();'
        return block

    code = FINALIZE_BLOCK_RE.sub(replace_finalize_block, code)

    # --- 3a) Remove input_file declaration
    code = INPUT_FILE_RE.sub('', code)

    # --- 3b + "all uses": remove ParmParse pp(...) with 4 args and any uses
    code = remove_parmparse_and_uses(code)

    # --- 4) setupAndRun(<string here>) -> setupAndRun()
    code = SETUP_AND_RUN_RE.sub('setupAndRun()', code)

    # --- 5) RefCountedPtr<Class> var = ... -> auto var = ...
    code = REFCOUNTEDPTR_LHS_RE.sub(r'auto \1 =', code)

    return code


def process_file(path: str) -> None:
    with open(path, 'r') as f:
        content = f.read()

    # Only bother with files that actually have main() — executables
    if 'main(' not in content:
        return

    new_content = transform_content(content)

    if new_content != content:
        with open(path, 'w') as f:
            f.write(new_content)
        print(f"Transformed: {path}")


def walk_cpp_files(root: str) -> None:
    for dirpath, dirnames, filenames in os.walk(root):
        for fname in filenames:
            if fname.endswith('.cpp'):
                process_file(os.path.join(dirpath, fname))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 refactor_chombo_execs_v2.py <directory>")
        sys.exit(1)
    walk_cpp_files(sys.argv[1])
