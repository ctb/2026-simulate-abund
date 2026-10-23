#! /usr/bin/env python
"""
Generate a mutated version of a FASTA file.
"""
import sys
import argparse
import random
import screed

DNA=["A", "C", "G", "T"]

def main():
    p = argparse.ArgumentParser()
    p.add_argument('fasta_file')
    p.add_argument('-m', '--mutation-p', type=int, default=4,
                   help="percentage of bases to change (default: 4)")
    p.add_argument('--new-prefix', required=False)
    p.add_argument('-o', '--output-fasta', required=True)
    p.add_argument('--seed', type=int, default=42)
    args = p.parse_args()

    random.seed(args.seed)

    records = [ record for record in screed.open(args.fasta_file) ]

    if args.new_prefix:
        for r in records:
            num = r.name.split('_')[1]
            new_name = f"{args.new_prefix}_{num}"
            r.name = new_name

    ofp = open(args.output_fasta, 'wt')

    for r in records:
        newseq = list(r.sequence)
        for i in range(int(round(len(newseq) * args.mutation_p / 100))):
            pos = random.randrange(0, len(newseq))
            oldch = newseq[pos]
            ch = oldch

            # pick a different base.
            while ch == oldch:
                ch = random.choice(DNA)

            newseq[pos] = ch.lower()

        newseq = "".join(newseq)
            
        ofp.write(f'>{r.name}\n{newseq}\n')

    ofp.close()


if __name__ == '__main__':
    sys.exit(main())

