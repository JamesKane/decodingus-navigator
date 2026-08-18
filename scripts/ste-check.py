"""Check comments against the mechanically-checkable ASD-STE100 Simplified Technical English rules.

Advisory. It always exits 0 -- it reports, it does not gate. See documents/STE-dictionary.md for
the project Technical Names and Technical Verbs, which this script cannot check.

  python3 scripts/ste-check.py                 # all Rust comments (the enforced scope)
  python3 scripts/ste-check.py --detail        # with examples
  python3 scripts/ste-check.py path/to/file.rs # one file
  python3 scripts/ste-check.py --all           # include Markdown (out of scope, for information)


Covers the rules a script can judge without a POS tagger:
  STE 1.x  approved vocabulary (curated subset of common offenders + replacements)
  STE 2.x  no -ing participles used as verbs/adjectives
  STE 3.x  active voice; simple tenses
  STE 5.x  noun clusters of more than three words
  STE 6.x  sentence length (20 procedural / 25 descriptive)
  STE 7.x  paragraph length (max 6 sentences)
  STE 8.x  no slang, idiom, metaphor, jargon
"""
import os, re, sys, json
from collections import Counter, defaultdict

# --- STE non-approved words -> approved alternative (curated, high-frequency subset) -------
NOT_APPROVED = {
    "utilize": "use", "utilise": "use", "utilizes": "uses", "utilized": "used",
    "via": "by / through", "per": "for each", "vice": "instead of",
    "prior to": "before", "subsequent to": "after", "in order to": "to",
    "due to": "because of", "owing to": "because of", "as well as": "and",
    "in the event that": "if", "in case": "if", "provided that": "if",
    "additionally": "also", "furthermore": "also", "moreover": "also",
    "however": "but", "nevertheless": "but", "nonetheless": "but",
    "hence": "so", "thus": "so", "therefore": "so", "whilst": "while",
    "amongst": "among", "whereby": "by which", "wherein": "in which",
    "aforementioned": "the ... described before", "said": "the",
    "obtain": "get", "acquire": "get", "commence": "start", "initiate": "start",
    "terminate": "stop", "cease": "stop", "endeavour": "try", "attempt": "try",
    "ascertain": "find out", "determine": "find", "require": "need",
    "sufficient": "enough", "numerous": "many", "multiple": "more than one",
    "approximately": "about", "regarding": "about", "concerning": "about",
    "possess": "have", "purchase": "buy", "assist": "help", "permit": "let",
    "indicate": "show", "demonstrate": "show", "illustrate": "show",
    "facilitate": "help", "leverage": "use", "handle": "control",
    "perform": "do", "conduct": "do", "execute": "do",
    "ensure": "make sure", "verify": "check", "validate": "check",
    "prohibit": "do not let", "eliminate": "remove", "modify": "change",
    "accomplish": "do", "encounter": "find", "identify": "find",
    "sole": "only", "solely": "only", "merely": "only", "simply": "only",
    "essentially": "", "basically": "", "actually": "", "really": "",
    "quite": "", "rather": "", "fairly": "", "somewhat": "",
    "various": "different", "several": "some", "certain": "some",
    "considerable": "large", "substantial": "large", "significant": "large",
    "minimal": "small", "optimal": "best", "optimum": "best",
    "prior": "earlier", "latter": "second", "former": "first",
    "cannot": "can not", "won't": "will not", "don't": "do not",
    "doesn't": "does not", "isn't": "is not", "it's": "it is",
    "we've": "we have", "they're": "they are", "that's": "that is",
    "wouldn't": "would not", "couldn't": "could not", "didn't": "did not",
    "hasn't": "has not", "haven't": "have not", "aren't": "are not",
    "let's": "let us", "there's": "there is", "what's": "what is",
}

# --- project Technical Names, read from documents/STE-dictionary.md -------------------------
def _technical_names():
    """Words the project declares as Technical Names / Technical Verbs. STE permits these."""
    names = set()
    try:
        here = os.path.dirname(os.path.abspath(globals().get("__file__", "scripts/x")))
        doc = open(os.path.join(here, os.pardir, "documents", "STE-dictionary.md"),
                   encoding="utf-8").read()
    except OSError:
        return names
    body = doc.split("## Technical Names", 1)
    if len(body) < 2:
        return names
    body = body[1].split("## Words to avoid", 1)[0]
    for tok in re.split(r"[\u00b7|\n,()]", body):
        tok = tok.strip().strip("*`").lower()
        if tok and re.fullmatch(r"[a-z][a-z0-9 -]*", tok):
            names.add(tok)
            for w in tok.split():
                names.add(w)
    return names


TECHNICAL = _technical_names()

# --- idiom / metaphor / informal (STE 8) ----------------------------------------------------
IDIOMS = [
    "God file", "god file", "under the hood", "out of the box", "rule of thumb",
    "hand-rolled", "hand rolled", "hammering", "hammer", "byte for byte",
    "byte-for-byte", "the headline", "falls out of", "fell behind", "pays off",
    "earns its keep", "earned its keep", "loose ends", "knock-on", "gotcha",
    "smell test", "cheap", "expensive", "dead code", "wired in", "baked in",
    "boils down", "kick off", "kicks off", "spin up", "tear down", "teardown",
    "in flight", "in the wild", "happy path", "sad path", "sanity check",
    "belt and braces", "first-class", "second-class", "nuked", "blew up",
    "silently", "quietly", "magic", "magical", "clever", "ugly", "nasty",
    "painful", "pain", "trivial", "obvious", "simply put", "of course",
    "arguably", "unfortunately", "sadly", "luckily", "surprisingly",
    "the point is", "worth noting", "note that", "keep in mind", "bear in mind",
    "a stone's throw", "on the fly", "at a glance", "hunting for", "chase",
    "drift apart", "drift", "stands in the way", "standing in for",
]

BE = r"(?:is|are|was|were|be|been|being|get|gets|got)"
# Irregular participles must match as whole words. As a suffix, "sent" matched "present" and
# "set" matched "offset", which reported plain adjectives as passive voice.
IRREGULAR = r"built|set|put|read|kept|held|left|made|sent|done|shown|known|thrown|written|driven"
PASSIVE = re.compile(rf"\b{BE}\s+(?:\w+ly\s+)?(\w+(?:ed|en|own|ung)|(?:{IRREGULAR}))\b", re.I)
# Words that end in -en/-ed but are not participles.
NOT_PARTICIPLE = {
    "often", "when", "then", "even", "open", "green", "seven", "ten", "children", "women", "men",
    "garden", "token", "golden", "sudden", "hidden", "wooden", "kitchen", "screen", "between",
    "need", "indeed", "speed", "seed", "feed", "exceed", "red", "bed", "led", "ahead", "instead",
}
ING = re.compile(r"\b(\w{3,}ing)\b", re.I)
ING_OK = {
    # technical names / gerund-nouns STE permits as established terms
    "string", "strings", "ring", "spring", "during", "bring", "thing", "things",
    "nothing", "something", "anything", "everything", "morning", "king", "wing",
    "sing", "swing", "engineering", "sequencing", "painting", "matching",
    "mapping", "encoding", "decoding", "indexing", "logging", "polling",
    "signing", "warning", "warnings", "setting", "settings", "heading",
    "listing", "ordering", "padding", "casing", "tracking", "caching",
    "processing", "pending", "missing", "remaining", "existing", "following",
    "corresponding", "underlying", "according", "including",
    # Technical Names from documents/STE-dictionary.md that end in -ing.
    "operating", "genotyping", "pacing", "sampling", "scaling", "streaming", "spilling", "phasing",
    "binning", "masking", "trimming", "clipping", "calling", "sorting", "merging",
    "reading", "writing", "counting", "timing", "build", "backing",
}
SENT_SPLIT = re.compile(r"(?<=[.!?:;])\s+(?=[A-Z`\[(])")


def sentences(text):
    return [s.strip() for s in SENT_SPLIT.split(text) if s.strip()]


def strip_code(text):
    """Remove inline code, fenced code, links and paths — STE judges prose, not identifiers."""
    text = re.sub(r"```.*?```", " ", text, flags=re.S)
    text = re.sub(r"`[^`]*`", " CODE ", text)
    text = re.sub(r"https?://\S+", " URL ", text)
    text = re.sub(r"\[\[?[^\]]*\]\]?(\([^)]*\))?", " LINK ", text)
    text = re.sub(r"\b[\w/.-]+\.(rs|md|sql|toml|json|yml)\b", " FILE ", text)
    return text


def extract_rust(path):
    """Comment text from a .rs file, as (line_no, text)."""
    out = []
    for i, raw in enumerate(open(path, encoding="utf-8", errors="replace"), 1):
        s = raw.strip()
        m = re.match(r"^(///|//!|//)\s?(.*)$", s)
        if m and not s.startswith("////"):
            out.append((i, m.group(2)))
    return out


def extract_md(path):
    out = []
    fence = False
    for i, raw in enumerate(open(path, encoding="utf-8", errors="replace"), 1):
        s = raw.rstrip("\n")
        if s.strip().startswith("```"):
            fence = not fence
            continue
        if fence or not s.strip() or s.strip().startswith(("|", ">", "#", "---")):
            continue
        out.append((i, s.strip()))
    return out


def analyse(items, kind):
    """items: [(line, text)] -> violations by rule."""
    v = defaultdict(list)
    # Group consecutive lines into paragraphs for sentence-level rules.
    para, start, last = [], None, None
    paras = []
    for ln, t in items:
        # A blank comment line, a gap in line numbers (code between doc blocks), or the start of a
        # Markdown list item ends a paragraph. A list item is its own unit of prose.
        if not t.strip() or (last is not None and ln != last + 1) or re.match(r"[-*+]\s|\d+\.\s", t.strip()):
            if para:
                paras.append((start, " ".join(para)))
                para, start = [], None
        if t.strip():
            if start is None:
                start = ln
            para.append(t)
        last = ln
    if para:
        paras.append((start, " ".join(para)))

    for ln, ptext in paras:
        clean = strip_code(ptext)
        sents = sentences(clean)
        if len(sents) > 6:
            v["STE7 paragraph >6 sentences"].append((ln, f"{len(sents)} sentences"))
        for s in sents:
            words = re.findall(r"[A-Za-z][\w'-]*", s)
            n = len(words)
            if n > 25:
                v["STE6 sentence >25 words"].append((ln, f"{n}w: {s[:90]}"))
            pm = PASSIVE.search(s)
            if pm and pm.group(1).lower() not in NOT_PARTICIPLE:
                v["STE3 passive voice"].append((ln, s[:90]))
            for w in ING.findall(s):
                if w.lower() not in ING_OK:
                    v["STE2 -ing form"].append((ln, w))
        low = " " + clean.lower() + " "
        for bad, good in NOT_APPROVED.items():
            if bad in TECHNICAL:
                continue
            if re.search(rf"\b{re.escape(bad)}\b", low):
                v["STE1 non-approved word"].append((ln, f"{bad} -> {good or 'delete'}"))
        for idiom in IDIOMS:
            if re.search(rf"\b{re.escape(idiom.lower())}\b", low):
                v["STE8 idiom/metaphor/informal"].append((ln, idiom))
        # Judge the code-stripped text: a fenced shell block puts cargo's `--` argument separator
        # in the paragraph, and that is not an em-dash aside.
        if "—" in clean or " -- " in clean:
            v["STE6 em-dash aside"].append((ln, ""))
    return v


def main():
    detail = "--detail" in sys.argv
    only = None
    for a in sys.argv[1:]:
        if not a.startswith("--"):
            only = a
    totals = Counter()
    per_file = Counter()
    examples = defaultdict(list)

    include_md = "--all" in sys.argv
    targets = []
    # Walk only the requested path, so a single-file or single-crate check is fast.
    root = only if only and os.path.isdir(only) else (os.path.dirname(only) or "." if only else ".")
    if only and os.path.isfile(only):
        root = os.path.dirname(only) or "."
    for dp, dn, fn in os.walk(root):
        if any(x in dp for x in ("/target", "/.git", "/node_modules", "/.claude/worktrees")):
            continue
        for f in fn:
            p = os.path.join(dp, f)
            if only and only not in p:
                continue
            if f.endswith(".rs"):
                targets.append((p, extract_rust, "rust"))
            elif f.endswith(".md") and include_md:
                targets.append((p, extract_md, "md"))

    for p, fn, kind in targets:
        try:
            v = analyse(fn(p), kind)
        except Exception:
            continue
        c = sum(len(x) for x in v.values())
        if c:
            per_file[p] = c
        for rule, hits in v.items():
            totals[rule] += len(hits)
            for h in hits[:2]:
                examples[rule].append((p, h))

    print(f"{'RULE':<34} {'VIOLATIONS':>10}")
    print("-" * 46)
    for rule, n in totals.most_common():
        print(f"{rule:<34} {n:>10,}")
    print("-" * 46)
    print(f"{'TOTAL':<34} {sum(totals.values()):>10,}")
    print(f"\nfiles with >=1 violation: {len(per_file)}")
    print("\nworst 15 files:")
    for p, n in per_file.most_common(15):
        print(f"  {n:>5}  {p}")
    if detail:
        print("\n=== examples ===")
        for rule in totals:
            print(f"\n## {rule}")
            for p, (ln, txt) in examples[rule][:6]:
                print(f"  {p}:{ln}  {txt}")


main()
