import defopt
import cyvcf2
from collections import defaultdict


def main(vcf_file: str, *, prefix: str = "error-counts"):
    """
    count errors (de novos) in a VCF file split by SNP/indel and by sample
    this expects particular annotations added for paper analyses

    :param vcf_file: path to VCF file
    :param prefix: prefix for output files
    """
    d = {'snp': defaultdict(lambda: defaultdict(int)),
         'indel': defaultdict(lambda: defaultdict(int))}

    lcr_d = {'snp': defaultdict(lambda: [0, 0]), 'indel': defaultdict(lambda: [0, 0])}

    counts_by_sample = defaultdict(int)

    vcf = cyvcf2.VCF(vcf_file)
    for variant in vcf:
        if variant.FILTER != 'PASS' and variant.FILTER is not None:
            continue

        error_samples = variant.INFO.get('dn').split(',')
        frag_count = variant.INFO.get('fraguracy_count', 0)
        lcr = bool(variant.INFO.get('LCR', 0))
        #frag_count = variant.INFO.get('fraguracy_samples', 0)
        snp = 'snp' if variant.is_snp else 'indel'
        for s in error_samples:
            d[snp][s][frag_count] += 1
            lcr_d[snp][s][lcr] += 1
            counts_by_sample[s] += 1

    max_count = 50
    # now print out the snp and indel dictionaries to separate files
    for snp in ['snp', 'indel']:
        with open(f'{prefix}.{snp}.errors', 'w') as f:
            f.write("sample\tsnp\tcutoff\tcount\n")
            for sample in d[snp]:
                # we want to count all the errors
                cur_keys = [c for c in d[snp][sample].keys() if c >= max_count]
                cum_sum = sum([d[snp][sample][c] for c in cur_keys])
                f.write(f'{sample}\t{snp}\t{max_count}\t{cum_sum}\n')
                for i in range(max_count - 1, -1, -1):
                    cum_sum += d[snp][sample].get(i, 0)
                    f.write(f'{sample}\t{snp}\t{i}\t{cumsum}\n')

        with open(f'{prefix}.{snp}.lcr', 'w') as f:
            f.write("sample\tsnp\tlcr\tcount\n")
            for sample in lcr_d[snp]:
                # we want to count all the errors
                f.write(f'{sample}\t{snp}\ttrue\t{lcr_d[snp][sample][True]}\n')
                f.write(f'{sample}\t{snp}\tfalse\t{lcr_d[snp][sample][False]}\n')
    print(dict(counts_by_sample))

if __name__ == '__main__':
    defopt.run(main)
