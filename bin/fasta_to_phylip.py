#!/usr/bin/env python

"""
Convert fasta alignments to relaxed phylip ones in constant memory.
Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this::
$ fasta_to_phylip seqs.fasta seqs.phylip
"""

import os
import random
import re
import string
import sys


class Sequence(object):
    """The Sequence object has a string *header* and
    various representations."""

    def __init__(self, header, seq):
        self.header = re.findall('^>(\S+)', header)[0]
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    @property
    def phylip(self):
        return self.header + " " + self.seq.replace('.','-') + "\n"

    @property
    def fasta(self):
        return ">" + self.header + "\n" + self.seq + "\n"


def fasta_parse(path):
    """Reads the file at *path* and yields
       Sequence objects in a lazy fashion"""
    header = ''
    seq = ''
    with open(path) as fasta_open:
        for line in fasta_open:
            line = line.strip('\n')
            if line.startswith('>'):
                if header: yield Sequence(header, seq)
                header = line
                seq = ''
                continue
            seq += line
    yield Sequence(header, seq)


def main():
    # Get the shell arguments #
    fasta_file = sys.argv[1]
    phylip_file = sys.argv[2]
    # Check that the path is valid #
    if not os.path.exists(fasta_file): raise Exception("No file at %s." % fasta_file)
    # Use our two functions #
    seqs = fasta_parse(fasta_file)
    # Write the output to temporary file #
    tm_path = phylip_file + '.' + ''.join(random.choice(string.ascii_letters) for i in range(10))
    # Count the sequences #
    count = 0
    with open(tm_path, 'w') as phylip_open:
        for seq in seqs:
            phylip_open.write(seq.phylip)
            count += 1
    # Add number of entries and length at the top #
    with open(tm_path, 'r') as old, open(phylip_file, 'w') as new:
        new.write(" " + str(count) + " " + str(len(seq)) + "\n")
        new.writelines(old)
    # Clean up #
    os.remove(tm_path)


if __name__ == "__main__":
    main()
