#!/usr/bin/env python3
"""
Apply static_cast fixes for bugprone-narrowing-conversions.
Read warnings from stdin (clang-tidy output) and apply static_cast.

Usage:
  ./run_clang_tidy.sh 2>&1 | python3 fix_narrowing.py
"""

import sys
import re
import os


def parse_warnings(lines):
    """Parse clang-tidy warning lines for narrowing conversions."""
    warnings = []
    for line in lines:
        m = re.match(
            r'^(.+?):(\d+):(\d+): warning: narrowing conversion from \'(.+?)\''
            r'.* to (?:signed type |)\'(.+?)\' .*\[bugprone-narrowing-conversions',
            line
        )
        if m:
            filepath, linenum, col, from_type, to_type = m.groups()
            if not os.path.isfile(filepath):
                continue
            warnings.append({
                'file': filepath,
                'line': int(linenum),
                'col': int(col),
                'from': from_type,
                'to': to_type,
            })
    return warnings


# Types that cannot be used standalone as a cast target
_UNRESOLVABLE_TYPES = {'mapped_type', 'value_type', 'size_type', 'difference_type',
                        'key_type', 'reference', 'pointer'}


def simplify_cast_type(to_type, from_type):
    """Map verbose type description to a clean cast type."""
    t = to_type.strip()
    # If the target type is a nested/associated container type we can't use directly,
    # try to infer from the source type
    if t in _UNRESOLVABLE_TYPES:
        # Check if source is size_t - target is probably int
        if 'unsigned' in from_type or 'size_t' in from_type:
            return 'int'
        return None  # Signal that we can't safely cast
    if re.search(r'\bReal\b', t) or ('double' in t and 'unsigned' not in t):
        return 'double'
    if 'float' in t:
        return 'float'
    if 'unsigned int' in t:
        return 'unsigned int'
    if 'unsigned long long' in t:
        return 'unsigned long long'
    if 'unsigned long' in t:
        return 'unsigned long'
    if re.search(r'\blong long\b', t):
        return 'long long'
    if 'long' in t:
        return 'long'
    if re.search(r'\bint\b', t):
        return 'int'
    return t


def find_expression_end(source, start):
    """
    Walk forward from `start` to find the end of a postfix-expression chain.

    Handles:
    - () [] for function calls and subscripts
    - <T> template arguments (when < follows identifier or >)
    - -> arrow operator (not a stop)
    - :: scope resolution (not a stop)

    Stops at depth 0 when encountering:
    - ; , ) ]
    - << >> (double angle = shift/stream operators)
    - ? : (ternary) — : only when not preceded by :
    - * / % (multiplicative) and + (additive, not unary)
    - Standalone < > (comparison) at depth 0 and not inside template brackets
    """
    depth_paren = 0
    depth_bracket = 0
    depth_template = 0
    in_string = False
    in_char = False
    escape = False

    def is_ident_char(c):
        return c.isalnum() or c == '_'

    i = start
    while i < len(source):
        c = source[i]
        nxt = source[i + 1] if i + 1 < len(source) else ''
        prev = source[i - 1] if i > 0 else ''

        if escape:
            escape = False
            i += 1
            continue
        if c == '\\' and (in_string or in_char):
            escape = True
            i += 1
            continue
        if c == '"' and not in_char:
            in_string = not in_string
            i += 1
            continue
        if c == "'" and not in_string:
            in_char = not in_char
            i += 1
            continue
        if in_string or in_char:
            i += 1
            continue

        # Track parentheses
        if c == '(':
            depth_paren += 1
        elif c == ')':
            if depth_paren == 0 and depth_template == 0:
                return i
            if depth_paren > 0:
                depth_paren -= 1

        # Track square brackets
        elif c == '[':
            depth_bracket += 1
        elif c == ']':
            if depth_bracket == 0:
                return i
            depth_bracket -= 1

        # Only process operators when inside no parens/brackets
        elif depth_paren == 0 and depth_bracket == 0:
            if c in (',', ';'):
                return i

            # Ternary ? always stops
            if c == '?':
                return i

            # : stops unless it's :: (scope resolution)
            if c == ':':
                if nxt != ':' and prev != ':':
                    return i

            # << (stream/shift) stops; single < starts template or stops
            if c == '<':
                if nxt == '<':
                    # << operator — stop
                    return i
                else:
                    # Could be template bracket: < after identifier or ) or >
                    if is_ident_char(prev) or prev in (')', ']', '>'):
                        depth_template += 1
                    else:
                        # Comparison < — stop
                        return i

            elif c == '>':
                if nxt == '>':
                    # >> operator — stop
                    return i
                elif depth_template > 0:
                    depth_template -= 1
                elif prev == '-':
                    pass  # part of -> arrow operator
                else:
                    # Comparison > — stop
                    return i

            elif c == '-':
                if nxt == '>' or nxt == '-':
                    pass  # -> arrow or -- decrement, part of expression
                else:
                    return i  # subtraction operator

            elif c in ('*', '/', '%'):
                return i

            elif c == '+':
                if nxt != '+':
                    return i  # addition (not ++)

        i += 1
    return i


def line_col_to_offset(source, line, col):
    """Convert 1-based line/col to 0-based byte offset."""
    offset = 0
    for i, ln in enumerate(source.split('\n'), 1):
        if i == line:
            return offset + col - 1
        offset += len(ln) + 1
    return offset


def apply_cast(source, line, col, cast_type):
    """Wrap expression at (line, col) with static_cast<cast_type>(...)."""
    pos = line_col_to_offset(source, line, col)
    if pos >= len(source):
        return None

    expr_start = pos
    expr_end = find_expression_end(source, pos)
    expr = source[expr_start:expr_end].rstrip()

    if not expr:
        return None
    if expr.startswith('static_cast'):
        return None  # Already wrapped

    replacement = f"static_cast<{cast_type}>({expr})"
    return source[:expr_start] + replacement + source[expr_end:]


def process_file(filepath, file_warnings):
    """Apply all narrowing cast fixes to a single file."""
    with open(filepath, 'r') as f:
        source = f.read()

    # Sort in reverse order so earlier fixes don't shift offsets for later ones
    file_warnings = sorted(file_warnings, key=lambda w: (w['line'], w['col']), reverse=True)

    changes = 0
    skipped = []
    for w in file_warnings:
        cast_type = simplify_cast_type(w['to'], w['from'])
        if cast_type is None:
            skipped.append(f"  SKIP {w['line']}:{w['col']} — unresolvable type ({w['from']} -> {w['to']})")
            continue
        new_source = apply_cast(source, w['line'], w['col'], cast_type)
        if new_source is not None and new_source != source:
            source = new_source
            changes += 1
        else:
            skipped.append(f"  SKIP {w['line']}:{w['col']} ({w['from']} -> {w['to']})")

    if changes > 0:
        with open(filepath, 'w') as f:
            f.write(source)

    return changes, skipped


if __name__ == '__main__':
    input_lines = sys.stdin.read().splitlines()
    warnings = parse_warnings(input_lines)
    print(f"Parsed {len(warnings)} narrowing-conversion warnings")

    # Group by file
    by_file = {}
    for w in warnings:
        by_file.setdefault(w['file'], []).append(w)

    total = 0
    all_skipped = []
    for filepath in sorted(by_file):
        changes, skipped = process_file(filepath, by_file[filepath])
        rel = os.path.relpath(filepath)
        print(f"  {rel}: {changes} fixes applied" + (f", {len(skipped)} skipped" if skipped else ""))
        all_skipped.extend([(filepath, s) for s in skipped])
        total += changes

    if all_skipped:
        print("\nSkipped locations (need manual review):")
        for fp, s in all_skipped:
            print(f"  {os.path.relpath(fp)}: {s}")

    print(f"\nTotal: {total} fixes applied")
