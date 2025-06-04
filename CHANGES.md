# v0.2.7
+ indels: report length and bq-bin
+ combine_counts: better error messages and handle NA

# v0.2.6
+ report cases where neither bases matches the reference as NN:1 when --reference-as-truth is passed.
+ add lua expressions to filter reads

# v0.2.5
+ respect include and exclude for indels and denominator calculation (thanks Jason for reporting)

# v0.2.4
+ Allow chromosomes longer than u8::MAX (#11 thanks @pontushojer for reporting)
+ Fix: when an include region was given and a non-seen chromosome was queried, it would return all intervals in that chromosome (#10 thanks Jason)

# v0.2.3

+ Add --chromosome option to restrict analysis to a single chromosome.
+ fix counts of indel errors
+ fix distance to homopolymer (really)


# v0.2.2

+ Fix distance to homopolymer
