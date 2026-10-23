#! /usr/bin/env python
"""
Make simulated (error free) reads from multiple genomes at some set of ratios.
"""
import sys
import argparse
import random
import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('fasta_files', nargs='+')
    p.add_argument('-C', '--coverage', type=float, default=100)
    p.add_argument('--mix', type=str)
    p.add_argument('-L', '--read-lengths', type=int, default=100)
    p.add_argument('-o', '--output-fasta', required=True)
    p.add_argument('--seed', type=int, default=42)
    args = p.parse_args()

    mixtures = [ int(x) for x in args.mix.split(':') ]
    total = sum(mixtures)
    fractions = [ x/total for x in mixtures ]
    print('mixtures:', fractions)
    assert len(fractions) == len(args.fasta_files)

    readlen = args.read_lengths
    random.seed(args.seed)

    total_wrote = 0
    ofp = open(args.output_fasta, 'wt')
    for frac, filename in zip(fractions, args.fasta_files):
        seqs = [ record.sequence for record in screed.open(filename) ]
        n_seqs = len(seqs)

        total = sum([ len(x) for x in seqs ])
        num_reads = int(round(total * args.coverage * frac / readlen))

        for n in range(num_reads):
            name = f"read{n}"
            seqi = random.randrange(0, n_seqs)
            seq = seqs[seqi]

            pos = random.randrange(0, len(seq) - readlen)
            read = seq[pos:pos+readlen]
            if random.randint(0, 1): # rc?
                read = screed.rc(read)

            ofp.write(f'>{name}\n{read}\n')
            total_wrote += len(read)

    ofp.close()
    print(f"wrote {total_wrote/1000:.1f} kb")


if __name__ == '__main__':
    sys.exit(main())

