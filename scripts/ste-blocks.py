"""Print each comment block that has an STE violation, whole and ready to rewrite.

`ste-check.py` counts violations and names the rule. This prints the *block* around each one, so a
rewrite has the full context and does not have to be reassembled from line numbers.

  python3 scripts/ste-blocks.py crates/navigator-app/src/lib.rs

The header of each block gives its line range and the rules that fired inside it. See
documents/STE-dictionary.md for the workflow this belongs to.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CHECK = os.path.join(HERE, "ste-check.py")

# ste-check.py runs its own main() at import, so exec the source with that call removed.
src = open(CHECK, encoding="utf-8").read().replace("\nmain()\n", "\n")
ns = {}
exec(compile(src, CHECK, "exec"), ns)

if len(sys.argv) < 2:
    sys.exit("usage: ste-blocks.py <file.rs>")
path = sys.argv[1]
lines = open(path, encoding="utf-8").read().splitlines()
violations = ns["analyse"](ns["extract_rust"](path), "rust")

# line -> the rules that fired on it, with a short excerpt where the rule carries one.
by_line = {}
for rule, hits in violations.items():
    for line_no, detail in hits:
        tag = rule.split()[0] + (f":{detail[:22]}" if detail else "")
        by_line.setdefault(line_no, set()).add(tag)


def is_comment(i):
    return 0 <= i < len(lines) and lines[i].strip().startswith(("///", "//!", "//"))


printed = set()
for line_no in sorted(by_line):
    first = last = line_no - 1
    while is_comment(first - 1):
        first -= 1
    while is_comment(last + 1):
        last += 1
    if (first, last) in printed:
        continue
    printed.add((first, last))
    tags = set()
    for i in range(first + 1, last + 2):
        tags |= by_line.get(i, set())
    print(f"===== {path}:{first + 1}-{last + 1}  [{', '.join(sorted(tags))}]")
    for i in range(first, last + 1):
        print(lines[i])
    print()

total = sum(len(h) for h in violations.values())
print(f"# {total} violations in {len(printed)} blocks")
