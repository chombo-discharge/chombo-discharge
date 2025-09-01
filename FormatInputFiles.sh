#!/bin/bash

# Thanks to ChatGPT for this one. Love ya. 

find . \( -name '*.inputs' -o -name '*.options' \) ! -path '*/Submodules/*' | while read -r file; do

  # Part 1: Align '=' characters across lines in each block
  awk '
  BEGIN {
    RS = "";
    FS = "\n";
  }
  {
    n = split($0, lines, "\n")

    max_left_len = 0
    for (i=1; i<=n; i++) {
      line = lines[i]
      if (match(line, /^[ \t]*#/)) continue

      eq_pos = -1
      for (j=1; j<=length(line); j++) {
        if (substr(line,j,1) == "=") {
          eq_pos = j
          break
        }
      }

      if (eq_pos > 0) {
        left_part = substr(line, 1, eq_pos-1)
        sub(/[ \t]+$/, "", left_part)

        left_len = length(left_part)
        if (left_len > max_left_len) max_left_len = left_len
      }
    }

    eq_align = max_left_len + 5

    for (i=1; i<=n; i++) {
      line = lines[i]

      if (match(line, /^[ \t]*#/)) {
        print line
        continue
      }

      eq_pos = -1
      for (j=1; j<=length(line); j++) {
        if (substr(line,j,1) == "=") {
          eq_pos = j
          break
        }
      }

      if (eq_pos > 0) {
        left_part = substr(line, 1, eq_pos-1)
        right_part = substr(line, eq_pos+1)

        sub(/[ \t]+$/, "", left_part)
        sub(/^[ \t]+/, "", right_part)

        pad_len = eq_align - length(left_part) - 1
        if (pad_len < 0) pad_len = 0

        left_part = left_part sprintf("%*s", pad_len, "")

        line = left_part " = " right_part
      }

      print line
    }
    print ""
  }
  ' "$file" > "$file.tmp" && mv "$file.tmp" "$file"


  # Part 2: Align "#" characters across lines in each block (include lines starting with # indented >5 spaces or tabs)
  awk '
  BEGIN {
    RS = "";
    FS = "\n";
  }
  {
    n = split($0, lines, "\n")

    max_left_text_len = 0
    for (i=1; i<=n; i++) {
      line = lines[i]

      if (match(line, /^[ \t]+#/)) {
        indent_len = 0
        tab_found = 0
        for (j=1; j<=length(line); j++) {
          c = substr(line, j, 1)
          if (c == " ") indent_len++
          else if (c == "\t") {
            tab_found = 1
            break
          } else {
            break
          }
        }
        if (indent_len <= 5 && !tab_found) continue
      } else if (match(line, /^[ \t]*#/)) {
        continue
      }

      hash_pos = index(line, "#")
      if (hash_pos > 0) {
        left_substr = substr(line, 1, hash_pos - 1)
        sub(/[ \t]+$/, "", left_substr)
        left_len = length(left_substr)
      } else {
        sub(/[ \t]+$/, "", line)
        left_len = length(line)
      }

      if (left_len > max_left_text_len) max_left_text_len = left_len
    }

    hash_align = max_left_text_len + 5

    for (i=1; i<=n; i++) {
      line = lines[i]

      if (match(line, /^[ \t]+#/)) {
        indent_len = 0
        tab_found = 0
        for (j=1; j<=length(line); j++) {
          c = substr(line, j, 1)
          if (c == " ") indent_len++
          else if (c == "\t") {
            tab_found = 1
            break
          } else {
            break
          }
        }
        if (indent_len <= 5 && !tab_found) {
          print line
          continue
        }
      } else if (match(line, /^[ \t]*#/)) {
        print line
        continue
      }

      hash_pos = index(line, "#")

      if (hash_pos > 0) {
        left_part = substr(line, 1, hash_pos-1)
        right_part = substr(line, hash_pos)

        sub(/[ \t]+$/, "", left_part)

        pad_len = hash_align - length(left_part) - 1
        if (pad_len < 0) pad_len = 0

        left_part = left_part sprintf("%*s", pad_len, "")

        line = left_part right_part
      }

      print line
    }
    print ""
  }
  ' "$file" > "$file.tmp" && mv "$file.tmp" "$file"


  # Part 3: Remove trailing whitespace at each line
  sed -i 's/[ \t]*$//' "$file"


  # Part 4: Replace all '##' with '#'
  awk '
  {
    line = $0

    while (index(line, "##") > 0) {
      sub("##", "#", line)
    }

    print line
  }
  ' "$file" > "$file.tmp" && mv "$file.tmp" "$file"


  # Part 5: Replace '#' by '##' in lines where '#' is present but not at line start
  awk '
  {
    line = $0
    if (match(line, /^[ \t]*#/)) {
      print line
    } else {
      if (index(line, "#") > 0) {
        sub("#", "##", line)
      }
      print line
    }
  }
  ' "$file" > "$file.tmp" && mv "$file.tmp" "$file"

  # Part 6: Remove empty line at the bottom of the file
  sed -i '${/^$/d;}' "$file"  

done
