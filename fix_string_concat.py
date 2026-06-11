#!/usr/bin/env python3
"""Transform chained string operator+ to .append() chain in throwParser* calls."""

import re
import sys


def split_top_level_plus(expr):
    """Split expr on + at top-level (not inside parens, brackets, or string literals)."""
    parts = []
    depth = 0
    in_string = False
    escape = False
    string_char = None
    start = 0
    i = 0
    while i < len(expr):
        c = expr[i]
        if escape:
            escape = False
            i += 1
            continue
        if c == '\\' and in_string:
            escape = True
            i += 1
            continue
        if in_string:
            if c == string_char:
                in_string = False
            i += 1
            continue
        if c in ('"', "'"):
            in_string = True
            string_char = c
            i += 1
            continue
        if c in ('(', '[', '{'):
            depth += 1
            i += 1
            continue
        if c in (')', ']', '}'):
            depth -= 1
            i += 1
            continue
        if c == '+' and depth == 0:
            parts.append(expr[start:i])
            start = i + 1
            i += 1
            continue
        i += 1
    parts.append(expr[start:])
    return parts


def strip_outer_whitespace_and_newlines(s):
    """Strip leading/trailing whitespace preserving internal structure."""
    return s.strip()


def transform_concat(arg_str):
    """Transform 'A + B + C' to 'std::string{A}.append(B).append(C)' if 2+ parts."""
    parts = split_top_level_plus(arg_str)
    # Strip whitespace from each part
    parts = [strip_outer_whitespace_and_newlines(p) for p in parts]
    # Filter empty parts (shouldn't happen with valid C++)
    parts = [p for p in parts if p]
    if len(parts) < 2:
        return None  # Nothing to transform
    # Build the .append() chain
    first = parts[0]
    rest = parts[1:]
    result = f"std::string{{{first}}}"
    for p in rest:
        result += f".append({p})"
    return result


def extract_call_arg(text, start_pos):
    """
    Given text and position of '(' after throwParserError/throwParserWarning,
    extract the full argument (until matching ')') and return (arg, end_pos).
    end_pos is the position just after the closing ')'.
    """
    assert text[start_pos] == '('
    depth = 0
    in_string = False
    escape = False
    string_char = None
    i = start_pos
    arg_start = start_pos + 1
    while i < len(text):
        c = text[i]
        if escape:
            escape = False
            i += 1
            continue
        if c == '\\' and in_string:
            escape = True
            i += 1
            continue
        if in_string:
            if c == string_char:
                in_string = False
            i += 1
            continue
        if c in ('"', "'"):
            in_string = True
            string_char = c
            i += 1
            continue
        if c == '(':
            depth += 1
            i += 1
            continue
        if c == ')':
            if depth == 1:
                return text[arg_start:i], i + 1
            depth -= 1
            i += 1
            continue
        i += 1
    return None, -1  # Unmatched paren


def transform_file(filename):
    with open(filename, 'r') as f:
        content = f.read()

    # Find all throwParserError( and throwParserWarning( occurrences
    pattern = re.compile(r'throwParser(?:Error|Warning)\s*\(')
    new_content = []
    pos = 0
    changes = 0

    for m in pattern.finditer(content):
        # Add content up to the match
        new_content.append(content[pos:m.start()])

        call_name = m.group(0)  # e.g. "throwParserError("
        paren_pos = m.end() - 1  # position of '('

        arg, end_pos = extract_call_arg(content, paren_pos)

        if arg is None:
            # Could not parse - emit original
            new_content.append(content[m.start():m.end()])
            pos = m.end()
            continue

        # Check if the arg has 2+ top-level + (3+ parts) — single + is fine
        parts_count = len(split_top_level_plus(arg))
        if parts_count >= 3:
            transformed = transform_concat(arg)
            if transformed is not None:
                # Reconstruct: call_name (without opening paren) + transformed + )
                func_part = call_name[:-1]  # remove trailing '('
                new_content.append(f"{func_part}({transformed})")
                pos = end_pos
                changes += 1
                continue

        # No transformation needed - emit original
        new_content.append(content[m.start():end_pos])
        pos = end_pos

    new_content.append(content[pos:])
    result = ''.join(new_content)

    if changes > 0:
        with open(filename, 'w') as f:
            f.write(result)
        print(f"{filename}: {changes} transformations applied")
    else:
        print(f"{filename}: no changes")

    return changes


if __name__ == '__main__':
    total = 0
    for fname in sys.argv[1:]:
        total += transform_file(fname)
    print(f"Total: {total} transformations")
