#! /usr/bin/env python
"""
Make simulated (error free) reads from a single genome.
"""
import sys
import argparse
import random
import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('fasta_file')
    p.add_argument('-C', '--coverage', type=float, default=100)
    p.add_argument('-L', '--read-lengths', type=int, default=100)
    p.add_argument('-o', '--output-fasta', required=True)
    p.add_argument('--seed', type=int, default=42)
    args = p.parse_args()

    random.seed(args.seed)

    seqs = [ record.sequence for record in screed.open(args.fasta_file) ]
    n_seqs = len(seqs)

    readlen = args.read_lengths
    total = sum([ len(x) for x in seqs ])
    num_reads = int(round(total * args.coverage / readlen))

    ofp = open(args.output_fasta, 'wt')

    total_wrote = 0
    for n in range(num_reads):
        name = f"read{n}"
        seqi = random.randrange(0, n_seqs)
        seq = seqs[seqi]

        # select a point in sequence at random
        pos = random.randrange(0, len(seq))

        # now, select an interval start from all possible interval starts
        # that include this point.
        # (this is called "rejection sampling".)
        lo = max(0, pos - readlen + 1)
        hi = min(pos, len(seq) - readlen)

        pos = random.randrange(lo, hi + 1)

        read = seq[pos:pos+readlen]
        if random.randint(0, 1): # rc?
            read = screed.rc(read)

        assert len(read) == readlen
        ofp.write(f'>{name}\n{read}\n')
        total_wrote += len(read)
    ofp.close()

    print(f"wrote {total_wrote/1000:.1f} kb")
    print(f"calculated coverage: {total_wrote / total:.1f}")
    print(f"desired coverage: {args.coverage:.1f}")


if __name__ == '__main__':
    sys.exit(main())

