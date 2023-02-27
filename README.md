# TODO

see locations like: samtools mpileup -r chr4:1801736-1801736 -q 10 -q 10 -d 100000 ./19610X1.mapped.sorted.bam
often have a lot of bases with 'g' clustered. Just ignore?

how to make pileup faster? do it and end when we have a list of variable sites?

output list of untrustworthy sites where there is an error in > 1 read.


# Fraguracy

`fraguracy` calculates error rates using overlapping reads in a fragment. This avoids some bias.
It does limit to the (potentially) small percentage of bases that overlap and it will sample less at the
beginning of read 1 and the end of read2 (which are known to have increased error rates).


## Bins

The aim is to create a model of errors. Many factors can be predictive of the likelihood of an error.
The dimensionality is a consideration because if the data is too sparse, prediction is less reliable.
Therefore we limit to: **Base-Quality**, **Mapping-Quality**, **Sequence Context**, **Read**, and **Position in Read**
as described and binned below. With those binnings we have 15,000 possible combinations (5 * 5 * 6 * 2 * 30 )

For each combination, while iterating over the bam, we store the number of errors and the number of total bases
in each bin. These become, respectively, the numerator and denominator for the error-rate for that set of parameters.

### Qualities (5)

Base-Qualities and Mapping Qualities will be binned to:

0. 0-5
1. 6-20
2. 21 - 40,
3. 41 - 59,
4. 60+

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



