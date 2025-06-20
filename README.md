<!--- 
# static binary build
RUSTFLAGS="-C target-feature=+crt-static" cargo build --release  --target x86_64-unknown-linux-gnu
--->
# Fraguracy

[![Rust](https://github.com/brentp/fraguracy/actions/workflows/rust.yml/badge.svg)](https://github.com/brentp/fraguracy/actions/workflows/rust.yml)

`fraguracy` calculates real error rates using overlapping paired-end reads in a fragment.
It reports a file of error positions and counts, along with a summary of errors by context, read-position, read-orientation (F or R) and base-quality.
While the overlap requirment does limit to the (potentially) small percentage of bases that overlap this can
still be useful to:

1. evaluate error rates within and among samples
2. find sites in the genome with high error rates
3. find data-driven cutoffs for allele fraction (`AF`) cutoffs in UMI or duplex sequencing.

# Usage

The `fraguracy` binary available in releases takes a bam or cram file and outputs error stats. The plotting is currently done via python.

```
$ fraguracy extract \
    --bin-size 1 \
    --output-prefix fraguracy-$sample- \
    --fasta $reference \
    $sample.bam [.. *.bam] \

$ python plot.py fraguracy-$sample-consensus-counts.txt # writes read.html

$ head fraguracy-$sample-errors.bed # records base position of every error observed and count of errors at that site.
chrom start stop bq_bin count contexts
chr1 75822283 75822284 05-19 6 AC:4,AT:2
chr1 75822287 75822288 20-36 4 TC:4
chr1 75822287 75822288 37-59 3 TC:3
chr1 75822287 75822288 60+ 2 CA:2
chr1 75822341 75822342 05-19 2 TC:2
chr1 75822352 75822353 20-36 2 GT:2
chr1 75822360 75822361 20-36 2 AG:2
chr1 241850751 241850752 37-59 2 TC:2
chr1 241850752 241850753 20-36 2 TA:1,TC:1
```

There is also an `$sample-indel-errors.bed` file that contains the columns:

```
chrom   start   stop    count
```

The errors files are useful to find **positions that are frequent errors** -- having count > 1 or with multiple bq_bins showing the same position.

If multiple samples are given (multiple bam files) then each sample is processed in parallel and $prefix-total-counts.txt and $prefix-total-errors.bed will
be created which sum all values for all samples.

The plot.py will create an interactive plot that looks like this:

![frag-plot](https://user-images.githubusercontent.com/1739/225074861-7b5098d1-b5e9-4bab-8971-0a278f182aaa.png)

**NOTE** that depending on the goal it can be useful to run `fraguracy extract` once, then exclude sites that are very frequent errors and re-run,
this will prevent a small percentage of sites (often around homopolymers) from dominating the error profile.

## CLI

```
error profile pair overlaps in bam/cram

Usage: fraguracy extract [OPTIONS] [BAMS]...

Arguments:
  [BAMS]...  

Options:
  -f, --fasta <FASTA>
          fasta for use with crams and/or to use as 'truth'
  -o, --output-prefix <OUTPUT_PREFIX>
          prefix for output files [default: fraguracy-]
  -C, --chromosome <CHROMOSOME>
          restrict analysis to this chromosome
  -r, --regions <REGIONS>
          restrict analysis to the regions given in this BED file
  -e, --exclude-regions <EXCLUDE_REGIONS>
          exclude from analysis the regions given in this BED file
  -l, --lua-expression <LUA_EXPRESSION>
          optional lua expression to filter reads. returns true to skip read. e.g. 'return read.flags.secondary or read.flags.supplementary'.
  -m, --max-read-length <MAX_READ_LENGTH>
          indicate the maximum read length in the alignment file [default: 151]
  -b, --bin-size <BIN_SIZE>
          parition the read into chunks/bins of this size [default: 3]
  -Q, --min-mapping-quality <MIN_MAPPING_QUALITY>
          only consider pairs where both reads have this mapping-quality or higher (good to leave this high) [default: 50]
  -c, --ci <CI>
          method for confidence interval calculation (see rust bpci crate) [default: agresti-coull] [possible values: agresti-coull, wald, wilson]
  -n, --no-denominator
          do not calculate denominator. This can shorten runtime.
  -H, --homopolymer-regex <HOMOPOLYMER_REGEX>
          regex for homopolymer sequence to consider if denominator is calculated[default: A{3,}|C{3,}|G{3,}|T{3,}] [default: A{3,}|C{3,}|G{3,}|T{3,}]
  -t, --reference-as-truth
          use reference base as 'truth'
  -h, --help
          Print help
```

### Lua Expressions

The `extract` sub-command allows lua expressions with `-l` that indicate whether to skip a read. See [lua-api.md](lua-api.md) for a full description
of how to use this.

### Combine

`fraguracy extract` can also be run per-sample and then errors can be combined with `fraguracy combine-errors`:

```
Usage: fraguracy combine-errors [OPTIONS] --fai-path <FAI_PATH> [ERRORS]...

Arguments:
  [ERRORS]...  path to error bed files from extract

Options:
  -f, --fai-path <FAI_PATH>        path for to fai (not fasta) file
  -o, --output-path <OUTPUT_PATH>  path for output bed file [default: fraguracy-combined-errors.bed]
  -h, --help                       Print help
```

The output is a single file with the error counts from each sample summed. And an additional column indicating the
number of samples that containing the error is reported.

> ⚠️ **Warning**: You must send either indel error files or snp error files, not both!

## Bins

The aim is to create a model of errors. Many factors can be predictive of the likelihood of an error.
The dimensionality is a consideration because if the data is too sparse, prediction is less reliable.
Because we determine accuracy by the mapping, it is best to require a high mapping-quality.
Therefore we limit to: **Base-Quality**, **Sequence Context**, **Read**, and **Position in Read**
as described and binned below. With those binnings we have **189,720** possible combinations (5 *6* 2 *2* $read_length / $bin-size * 31 )

For each combination, while iterating over the bam, we store the number of errors and the number of total bases
in each bin. These become, respectively, the numerator and denominator for the error-rate for that set of parameters.

### Qualities (5)

Base-Qualities and Mapping Qualities will be binned to:

0. 0-5
1. 6-19
2. 20 - 36,
3. 37 - 59,
4. 60+

This means that the quantized base-qualities from nova-seq (2, 12, 23 and 37) are each in separate bins.
And other base-quality schemes are also paritioned sanely.

### Sequence Context (6)

0. C->A (G->T)
1. C->G (G->C)
2. C->T (G->A)
3. T->A (A->T)
4. T->C (A->G)
5. T->G (A->C)

### Read (2)

Read 1 or Read 2

### Read Position (50)

read position is simply divided by 3. so bins of 3 bases.

### Homopolymer distance (30)

The errors are also partitioned by homopolymer distance up to +- 15. all errors beyond 15 are put in the 15 base bin

# vcfanno

To use the errors files with vcfanno:

```
bgzip fraguracy/fraguracy-19610X19-errors.bed
tabix fraguracy/fraguracy-19610X19-errors.bed.gz

echo '
[[annotation]]
file="fraguracy/fraguracy-19610X19-errors.bed.gz"
columns=[4, 5]
names=["frag_bq_bin", "frag_errors"]
ops=["first", "first"]
' > conf.toml

vcfanno conf.toml $vcf > annotated.vcf # annotated.vcf will have entries for `frag_bq_bin` and `frag_errors` where there was an error found that was also a variant in the VCF.
```

## indel errors

A command like: `fraguracy extract -f $fasta -o $prefix $bam` will create the files needed to evaluate the indel error rate. To plot it, then use:

```
python scripts/analyze_indel_errors.py ${prefix}-indel-errors.bed.gz ${prefix}-counts.txt
```

Which will make a plot like this one:

![Image](https://github.com/user-attachments/assets/35c10634-a54b-48a2-9c5f-2c6363def26f)
