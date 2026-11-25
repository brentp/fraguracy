# Lua API Reference

This document describes the Lua API for filtering BAM records and overlapping bases in `fraguracy extract`.

There are two Lua filter flags:

- `-l` / `--lua-expression`: **Read-level** filter evaluated once per overlapping read pair (against read `a` and read `b`). Has access to `read`. Returning `true` skips the entire read pair.
- `-L` / `--lua-base-expression`: **Per-base** filter evaluated only for overlapping pairs that pass `-l`. Has access to both `read` and `base`. For overlapping base positions, returning `true` for either mate skips that position. For indel errors, it is evaluated only on the mate carrying the indel (`base.is_indel == true`). Keeping read-level checks in `-l` and position-specific checks in `-L` lets you opt into per-base evaluation only when needed.

## Flags

The `flags` object represents BAM flags returned from `read.flags` and provides the following fields:

- `paired`: Returns true if the read is paired
- `proper_pair`: Returns true if the read is in a proper pair
- `unmapped`: Returns true if the read is unmapped
- `mate_unmapped`: Returns true if the mate is unmapped
- `reverse`: Returns true if the read is on the reverse strand
- `forward`: Returns true if the read is on the forward strand
- `mate_reverse`: Returns true if the mate is on the reverse strand
- `mate_forward`: Returns true if the mate is on the forward strand
- `read_1`: Returns true if this is read 1
- `read_2`: Returns true if this is read 2
- `secondary`: Returns true if this is a secondary alignment
- `primary`: Returns true if this is the primary alignment
- `qcfail`: Returns true if the read fails quality checks
- `duplicate`: Returns true if the read is a duplicate
- `supplementary`: Returns true if this is a supplementary alignment
- `flag`: Returns the raw integer flag value

## Read (`-l` and `-L`)

The `read` object provides access to BAM record data with the following fields and methods:

### Fields

All of these are properties on the `read`, e.g. `read.mapping_quality`

- `mapping_quality`: Returns the mapping quality
- `flags`: Returns a `Flags` object
- `tid`: Returns the reference sequence ID
- `start`: Returns the 0-based start position
- `stop`: Returns the end position based on CIGAR
- `length`: Returns the sequence length
- `insert_size`: Returns the insert size
- `qname`: Returns the query name as a string
- `sequence`: Returns the read sequence as a string
- `soft_clips_3_prime`: Returns the number of soft-clipped bases at the 3' end
- `soft_clips_5_prime`: Returns the number of soft-clipped bases at the 5' end
- `base_counts`: Returns a table with counts of A, C, G, T, N in the read
- `n_proportion`: Returns the proportion of N bases in the read (see methods to limit to 3' or 5')
- `indel_count`: Returns the number of indels in the read
- `average_base_quality`: Returns the average base-quality in a read

### Methods

- `read:tag(tag_name)`: Returns the value of the specified BAM tag
- `read:n_proportion_3_prime(n:number)`: Returns the proportion of N bases within `n` of the 3' end of the read
- `read:n_proportion_5_prime(n:number)`: Returns the proportion of N bases within `n` of the 5' end of the read

## Base (`-L` / `--lua-base-expression` only)

The `base` object provides fast per-position access for each overlapping base (and indel error):

- `base.qpos`: 0-based query position in the read (for deletions, the query position of the base immediately following the deletion)
- `base.gpos`: 0-based genomic reference position
- `base.bq`: Phred base quality at `base.qpos` (for indels, the reported indel quality: the lower of the indel-carrying mate's quality at `base.qpos` and the next indel on the other mate if one exists)
- `base.distance_from_5prime`: 0-based distance from the 5' end of the read (`0` at 5' end, `length - 1` at 3' end)
- `base.distance_from_3prime`: 0-based distance from the 3' end of the read (`0` at 3' end, `length - 1` at 5' end)
- `base.base`: Base call at `base.qpos` (`"A"`, `"C"`, `"G"`, `"T"`, `"N"`)
- `base.is_indel`: `true` when evaluating an indel error (evaluated on the mate carrying the indel), `false` for base positions
- `base.indel_length`: `0` for base positions, positive for insertions, negative for deletions

> **Performance tip**: Put read-level predicates in `-l` and only position-dependent predicates in `-L`. When `-L` only accesses `base.*` (not `read.*` or `base.base`), fraguracy evaluates an entire overlap chunk inside Luau with zero per-base Rust callbacks.

## Usage Examples

```bash
# Read-level only (-l): skip pairs where either read has MQ < 20 or >10% N in last 10 bases
fraguracy extract -f ref.fa \
  -l 'return read.mapping_quality < 20 or read:n_proportion_3_prime(10) > 0.1' \
  sample.bam

# Combine read-level (-l) and per-base (-L): skip reads with avg BQ < 20, and skip bases within 5 bp of 3' end or with BQ < 15
fraguracy extract -f ref.fa \
  -l 'return read.average_base_quality < 20' \
  -L 'return base.distance_from_3prime < 5 or base.bq < 15' \
  sample.bam
```
