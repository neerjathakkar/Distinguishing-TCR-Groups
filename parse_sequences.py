# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

# Method for parsing a list of cdr sequences

import itertools


def txttoseqlist(filename, truncate_reads=False, num_reads=None, trim_first=True, rawtext=False, remove_special=False):
    # read in sequences
    with open(filename) as f:
        if rawtext:
            seq_list = f.readlines()
            seqs = ([x.splitlines() for x in seq_list])
            seqs = list(itertools.chain.from_iterable(seqs))
        else:
            seqs = f.readlines()
        seqs = [x.strip() for x in seqs]

    if not num_reads:
        num_reads = len(seqs)

    # get num_reads cdrs
    if truncate_reads:
        seqs = seqs[:num_reads]

    for i in range(len(seqs)):
        if trim_first:
            # remove the "C" from the beginning of each cdr
            seqs[i] = seqs[i][1:]

        # change all ~ to *
        seqs[i] = seqs[i].replace("~", "*")

        if remove_special:
            seqs[i]= seqs[i].replace("*", "")

    # remove duplicate sequences
    seqs = list(set(seqs))

    return seqs

