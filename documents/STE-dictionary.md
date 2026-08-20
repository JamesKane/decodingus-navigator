# Simplified Technical English — project dictionary

Navigator's comments and documentation follow **ASD-STE100 Simplified Technical English**.

STE permits three classes of word: the ~900 words in the STE approved dictionary, **Technical
Names**, and **Technical Verbs**. A project must declare its own Technical Names and Technical
Verbs, because the approved dictionary contains no domain vocabulary. This file is that
declaration. A word that is not in the STE dictionary and not in a table below must not appear in a
comment or a document.

`scripts/ste-check.py` enforces the mechanical rules. It cannot enforce vocabulary beyond the
curated list it carries, so use this file when you write.

## The rules, in short

| Rule | Requirement |
|---|---|
| STE 1 | Use one approved word for one meaning. Do not use a word as more than one part of speech. |
| STE 2 | Do not use an `-ing` form as a verb or an adjective. A Technical Name that ends in `-ing` is permitted. |
| STE 3 | Write in the active voice. Use the simple present, the simple past, or the simple future. |
| STE 4 | Write instructions as commands. Give one instruction in one sentence. |
| STE 5 | Do not use more than three nouns together. |
| STE 6 | Write sentences of 20 words or fewer in a procedure, and 25 words or fewer in a description. |
| STE 7 | Write paragraphs of six sentences or fewer. Write about one topic in one paragraph. |
| STE 8 | Do not use slang, idiom, metaphor, or humour. |
| STE 9 | Keep the articles (`a`, `an`, `the`). Do not remove words to make a sentence short. |

## Technical Names

A Technical Name is a noun. It can be a compound noun. It cannot be a verb.

### Genetics and sequence data

alignment · allele · ancestry · admixture · autosome · base · biosample · build · call · caller ·
chromosome · consensus · contig · coverage · depth · donor · genome · genotype ·
genotyping array · haplogroup · haplotype · indel · kit · lineage · marker · panel · pedigree ·
ploidy · position · read ·
read metrics · reference · reference genome · region · segment · sequence · sequence run · sex ·
signature · site · subject · variant · Y-STR

### File and data formats

BAM · BED · CRAM · FASTA · gVCF · index · JSON · masterVar · sidecar · TSV · VCF

### Application concepts

app · artifact · cache · command · event · liftover · migration · outbox · profile · project ·
query · realignment · record · row · schema · store · table · workspace · worker

### Local LLM

chat completion · context · fact sheet · model · model server · narration · prompt ·
reasoning model · token

### Federation

AppView · attestation · consent · device key · DID · exchange · handle · IBD · PDS · record key ·
session · signature · suggestion · token

### Technical Names that end in `-ing`

STE 2 forbids an `-ing` form as a verb or an adjective. These are declared nouns, so they are
permitted: genotyping array · mapping · sequencing · painting · matching · encoding · decoding ·
indexing · logging · polling · signing · setting · heading · listing · ordering · padding · casing ·
tracking · caching · processing · operating system · pacing · sampling · scaling · streaming ·
spilling · phasing · binning · masking · trimming · clipping · calling · sorting · merging ·
reasoning model · streaming · copying

## Technical Verbs

STE permits a project to declare Technical Verbs when no approved verb has the meaning. Use these
only with the meaning given.

| Verb | Meaning in this project |
|---|---|
| to align | to map reads to a reference genome |
| to cache | to keep a result for later use |
| to call | to find a genotype or a haplogroup from read data |
| to genotype | to find the alleles at a site |
| to import | to read a file into the workspace |
| to index | to make an index file for an alignment |
| to lift over | to change coordinates from one build to another |
| to publish | to send a record to a PDS |
| to realign | to map the reads of an alignment to a different reference genome |
| to sign | to make a cryptographic signature |
| to sync | to send records to a PDS and to read records from a PDS |

## Words to avoid, and what to write

The approved dictionary gives one word for one meaning. These replacements occur most frequently in
this codebase.

| Do not write | Write |
|---|---|
| additionally, furthermore, moreover | also |
| approximately | about |
| cannot | can not |
| due to, owing to | because of |
| ensure | make sure |
| however, nevertheless | but |
| in order to | to |
| indicate, demonstrate | show |
| multiple, numerous | many, more than one |
| perform, execute, conduct | do |
| prior to | before |
| require | need |
| simply, merely, solely | only, or delete |
| subsequent to | after |
| sufficient | enough |
| thus, hence, therefore | so |
| utilize, leverage | use |
| various | different |
| verify, validate | check |
| via | by, through |
| whilst | while |

Do not use a contraction. Write `do not`, not `don't`.

## What STE does not change

- **Identifiers.** Function, type, field, and file names are code. STE does not apply to them.
- **Code in a comment.** Text in backticks is code, not prose.
- **Commit messages and pull request text.** These are a record of a decision, not product
  documentation. They are outside the scope of this standard.
- **The rationale itself.** STE controls *how* you write a reason. It does not tell you to delete
  the reason. A comment must still say why the code is as it is.

## How to convert a file

This is the method that converted `navigator-resource` and 31 of the 33 files of `navigator-app`.
Follow it, and a file needs about three passes.

### The two tools

```bash
python3 scripts/ste-check.py  crates/navigator-app/src/lib.rs   # the count, by rule
python3 scripts/ste-blocks.py crates/navigator-app/src/lib.rs   # each block, whole
```

`ste-check.py` gives the score. `ste-blocks.py` gives the text to rewrite, with the rules that fired
on each block. Both accept a file, a directory, or no argument, which reads each Rust file. A single
file returns in about 0.03 seconds, so run the check after each edit.

### The loop

1. Run `ste-blocks.py` on the file and read the first 10 to 15 blocks.
2. Rewrite them in one batch. Use a Python script with a list of `(old, new)` pairs and
   `str.replace`. Do not edit by hand, and do not use a regular expression on prose. An exact string
   pair either matches or reports a miss, and it can not damage a line that you did not read.
3. Run `ste-check.py` on the file.
4. Repeat until the count is 0.

Write a small guard into the batch script, so a pair that no longer matches is loud:

```python
for a, b in subs:
    if a not in s:
        print("MISS", a[:55].replace("\n", " "))
        continue
    s = s.replace(a, b)
```

### Expect three passes, not one

The count does not go to zero in one pass, and that is normal. A first pass on a file with 50
violations usually leaves 15 to 25.

The reason is mechanical: you replace one 40-word sentence with three explanatory sentences, and one
of those three is 26 words. Rule 6 then fires on a line that did not exist before. The process
converges, and each pass is smaller than the last.

### What a rewrite must keep

Keep every fact and every reason. The measured numbers, the dates, the sample ids, the failure that
caused the code — each is the value of the comment. `blocktree.rs` in `navigator-app` is the worked
example: almost every constant there is a limit that a measurement defends, and each defence
survived the conversion.

Lose the compression and the voice. That is the trade this standard makes, and it is intended.

### Common false positives

The checker is careful but not perfect. Three cases came up often:

- **A domain term reads as an error.** `handle`, `panel`, `build`, and `genotyping array` are
  Technical Names. Add the term to this file, and the checker reads it from here. Do not rewrite
  the sentence around a term that is correct.
- **An `-ing` Technical Name.** `mapping`, `sequencing`, `reasoning model`. Add it to the `-ing`
  section above and to `ING_OK` in the checker.
- **A word inside a longer word.** The passive-voice rule once matched `present` as `pre` + `sent`.
  If a violation makes no sense, read the checker before you rewrite the prose.

### Verify

`cargo check -p <crate> --all-targets` after each batch of files. A doc comment can break the build:
a `[link]` that no longer resolves, or a line that starts with `-` or `+`, which
`clippy::doc_lazy_continuation` reads as a Markdown list.

**Separate a paragraph with a `///` line, never with a blank line.** A blank line ends the doc
comment. The text above it then documents nothing, which is error E0585 and stops the build. This is
easy to do when you split one long paragraph into three, and some files hold two doc comments that
ran together into one block, which invites the mistake. Scan for it after each batch:

```python
import glob
for f in glob.glob("crates/**/*.rs", recursive=True):
    L = open(f, encoding="utf-8").read().split("\n")
    for i in range(len(L) - 2):
        if L[i].strip().startswith("///") and not L[i + 1].strip() and L[i + 2].strip().startswith("///"):
            print(f"{f}:{i + 1}")
```

`cargo fmt --all` before each commit. The pre-commit hook enforces it.

### Two doc comments that ran together

Some doc blocks hold the documentation of **two different items**. One function's doc sits
physically on the next function, and the function it describes has none. The conversion finds these,
because a block that changes topic in the middle is easy to see when you rewrite it.

Do not repair one by hand while you convert. Note it, finish the file, and move every stranded doc
in one commit of its own. A doc comment that moves to a different item changes the API
documentation, and that is a different review from a change to the English.

This query finds the candidates. It reports each multi-paragraph doc that a documented item carries,
where the **next** item at the same indent has no doc at all:

```python
import re, glob
ITEM = re.compile(r'^(\s*)(?:pub(?:\([^)]*\))?\s+)?(?:async\s+)?(?:unsafe\s+)?(?:const\s+)?'
                  r'(?:fn|struct|enum|type|trait|static|const)\s+([A-Za-z_][A-Za-z0-9_]*)')
```

Read the first paragraph of each hit against the two item names. Most hits are a documented item
followed by a private helper, and those are correct as they stand. A hit is real only when the first
paragraph describes the *next* item. Of 64 candidates in this repo, 11 were real.

Three shapes turn up, and each needs a different repair:

- **A stranded doc.** The text describes a real item that has no doc. Move it there.
- **Two summaries of one item.** Somebody rewrote a doc and kept both versions. Merge them.
- **A doc for code that no longer exists.** Delete it. `haplogroup.rs` carried the doc of a
  per-alignment ancestry estimator that the consensus path replaced. It documented nothing, and it
  made the function below it read as though it did two jobs.

The relocation must not lose a word. `git diff --numstat` proves that: for a pure move the insertion
and deletion counts are equal, and every other difference must have a reason you can name.

### Where the work stands

Converted to zero: `navigator-resource`, `navigator-store/src/sig_cache.rs`, and 31 of the 33 files
of `navigator-app`.

Remaining in `navigator-app`: `src/lib.rs` (the type documentation at the top is converted, the
`impl App` body is not) and `src/haplogroup.rs`.

Not started: `navigator-analysis`, `navigator-ui`, `navigator-domain`, `navigator-align`,
`navigator-panelbuild`, `navigator-store` beyond `sig_cache`, `navigator-refgenome`,
`navigator-sync`. Run `python3 scripts/ste-check.py` for the current count.
