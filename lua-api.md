# Lua API Reference

This document describes the Lua API for filtering BAM records.

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

## Read

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

### Not Implemented

The following are not implemented as they require a per-base approach
which is not used in fraguracy.

- `bq` [NO]: Returns the base quality at the current position
- `distance_from_5prime` [NO]: Returns the distance from the 5' end of the read
- `distance_from_3prime` [NO]: Returns the distance from the 3' end of the read
- `qpos()` [NO]: Returns the query position

## Usage Example

```lua
-- skip reads with mapping quality >= 20 and not supplementary and where the proportion of N's in the last 10 bases is > 0.1
return read.mapping_quality >= 20 and not read.flags.supplementary and read:n_proportion_3_prime(10) > 0.1
```
