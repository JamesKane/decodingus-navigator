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
chromosome · consensus · contig · coverage · depth · donor · genome · genotype · haplogroup ·
haplotype · indel · kit · lineage · marker · panel · pedigree · ploidy · position · read ·
read metrics · reference · reference genome · region · segment · sequence · sequence run · sex ·
signature · site · subject · variant · Y-STR

### File and data formats

BAM · BED · CRAM · FASTA · gVCF · index · JSON · masterVar · sidecar · TSV · VCF

### Application concepts

app · artifact · cache · command · event · liftover · migration · outbox · profile · project ·
query · realignment · record · row · schema · store · table · workspace · worker

### Federation

AppView · attestation · consent · device key · DID · exchange · handle · IBD · PDS · record key ·
session · signature · suggestion · token

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
