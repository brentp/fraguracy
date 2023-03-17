# Fraguracy

`fraguracy` calculates real error rates using overlapping paired-end reads in a fragment. This avoids some bias.
It does limit to the (potentially) small percentage of bases that overlap and it will sample less at the
beginning of read 1 and the end of read2 (which are known to have increased error rates).

# Usage

The `fraguracy` binary available in releases takes a bam file (cram supported soon) and outputs error stats. The plotting is currently done via python.

```
$ fraguracy extract \
    --bin-size 1 \
    --output-prefix fraguracy-$sample-consensus- \
    --fasta $reference \
    $sample.bam [.. *.bam] \

$ python plot.py fraguracy-$sample-consensus-counts.txt # writes read.html

$ head fraguracy-$sample-errors.bed # records base position of every error observed and count of errors at that site.
chrom	start	stop	bq_bin	count
chr1	75822283	75822284	05-19	6	
chr1	75822287	75822288	20-36	4	
chr1	75822287	75822288	37-59	3	
chr1	75822287	75822288	60+	2	
chr1	75822341	75822342	05-19	2	
chr1	75822352	75822353	20-36	2	
chr1	75822360	75822361	20-36	2	
chr1	241850751	241850752	37-59	2	
chr1	241850752	241850753	20-36	2	
```

The errors file is useful to find **positions that are frequent errors** -- having count > 1 or with multiple bq_bins showing the same position.

If multiple samples are given (multiple bam files) then each sample is processed in parallel and $prefix-total-counts.txt and $prefix-total-errors.bed will
be created which sum all values for all samples.

The plot.py will create an interactive plot that looks like this:

![frag-plot](https://user-images.githubusercontent.com/1739/225074861-7b5098d1-b5e9-4bab-8971-0a278f182aaa.png)


## Bins

The aim is to create a model of errors. Many factors can be predictive of the likelihood of an error.
The dimensionality is a consideration because if the data is too sparse, prediction is less reliable.
Because we determine accuracy by the mapping, it is best to require a high mapping-quality.
Therefore we limit to: **Base-Quality**, **Sequence Context**, **Read**, and **Position in Read**
as described and binned below. With those binnings we have 15,000 possible combinations (5 * 6 * 2 * $read_length / $bin-size )

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

