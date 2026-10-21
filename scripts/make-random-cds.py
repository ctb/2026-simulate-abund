#! /usr/bin/env python
import sys
import argparse
import random


DNA=["A", "C", "G", "T"]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('prefix')
    p.add_argument('-N', '--num-sequences', type=int, default=100)
    p.add_argument('-L', '--sequence-lengths', type=int, default=1000)
    p.add_argument('-o', '--output-fasta', required=True)
    p.add_argument('--seed', type=int, default=42)
    args = p.parse_args()

    random.seed(args.seed)

    ofp = open(args.output_fasta, 'wt')

    for n in range(args.num_sequences):
        name = f"{args.prefix}_{n+1}"
        seq = DNA * (args.sequence_lengths // 4)
        random.shuffle(seq)
        seq = "".join(seq)

        ofp.write(f'>{name}\n{seq}\n')

    ofp.close()


if __name__ == '__main__':
    sys.exit(main())

