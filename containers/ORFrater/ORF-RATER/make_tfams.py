#! /usr/bin/env python

import argparse
from plastid.readers.bed import BED_Reader
from plastid.genomics.roitools import SegmentChain, positionlist_to_segments
from collections import defaultdict
import os
import sys
from time import strftime

parser = argparse.ArgumentParser(description='Identify overlapping transcripts from a bed file, termed transcript families ("tfams"). Each tfam will '
                                             'be assigned a name, either based on the transcript IDs or from an optional gene name file. It is '
                                             'recommended that this script be run in a newly created folder containing only the input BED file, '
                                             'generated by prune_transcripts.py, and possibly the GENENAMES file. It is also recommended that the '
                                             'default value ("tfams") be used for parameter TFAMSTEM, for consistency with later scripts.')

parser.add_argument('-g', '--genenames', help='Path to file with gene name lookup table. Formatted as a simple two-column tab-delimited file, with '
                                              'one transcript ID followed by the corresponding gene name on each line. Gene names may be repeated, '
                                              'though assigning the same name to non-overlapping transcripts will trigger a warning. Not every '
                                              'transcript must be assigned a gene name. If no file is provided, or if no gene names are available '
                                              'for any of the transcripts in a family, transcript IDs will be used as names.')
parser.add_argument('--inbed', default='transcripts.bed', help='Transcriptome BED-file (Default: transcripts.bed)')
parser.add_argument('--tfamstem', default='tfams', help='Output filestem. OUTSTEM.txt will be a tab-delimited file indicating which transcripts are '
                                                        'in which tfam. OUTSTEM.bed will be a bed file showing the genomic positions of each tfam. '
                                                        '(Default: tfams)')
parser.add_argument('-v', '--verbose', action='store_true', help='Output a log of progress and timing (to stdout)')
parser.add_argument('-f', '--force', action='store_true', help='Force file overwrite')
opts = parser.parse_args()

outbedname = opts.tfamstem + '.bed'
outtxtname = opts.tfamstem + '.txt'
if not opts.force:
    if os.path.exists(outbedname):
        raise IOError('%s exists; use --force to overwrite' % outbedname)
    if os.path.exists(outtxtname):
        raise IOError('%s exists; use --force to overwrite' % outtxtname)

if opts.verbose:
    sys.stdout.write(' '.join(sys.argv) + '\n')

    def logprint(nextstr):
        sys.stdout.write('[%s] %s\n' % (strftime('%Y-%m-%d %H:%M:%S'), nextstr))
        sys.stdout.flush()

    logprint('Identifying transcript overlaps')

tfams = {}  # will contain integer keys to a tuple: ([list of tid],(chrom,strand),{set of gcoord})
genlookup = defaultdict(dict)  # indexed by (chrom,strand) keys to an integer key to tfam
next_tfam = 0
processed = 0
with open(opts.inbed) as inbed:
    for trans in BED_Reader(inbed):
        pos_set = trans.get_position_set()
        currfams = {genlookup[(trans.chrom, trans.strand)][pos] for pos in pos_set if pos in genlookup[(trans.chrom, trans.strand)]}
        if currfams:
            newfam = min(currfams)
            currfams.discard(newfam)
            for fam in currfams:
                oldfam_item = tfams.pop(fam)
                tfams[newfam][0].extend(oldfam_item[0])  # extend the list of tids
                assert (tfams[newfam][1] == oldfam_item[1])  # chrom,strand should match
                tfams[newfam][2].update(oldfam_item[2])  # union the gcoords
            tfams[newfam][0].append(trans.attr['ID'])
            tfams[newfam][2].update(pos_set)
        else:
            newfam = next_tfam
            next_tfam += 1
            tfams[newfam] = ([trans.attr['ID']], (trans.chrom, trans.strand), pos_set)
        for pos in tfams[newfam][2]:
            genlookup[(trans.chrom, trans.strand)][pos] = newfam  # override old families and expand as needed
        processed += 1


def _choose_name(names):
    """Somewhat silly function that chooses a gene name if more than one are given for this gene. Chooses the shortest if that's unique, then one
    consisting only of letters and numbers if that's unique, then the one with fewest numbers in it, and finally the first by alphabetical order"""
    remaining = list(set(names))
    if len(remaining) > 1:
        minlen = len(remaining[0])
        nextset = [remaining[0]]
        for name in remaining[1:]:
            if len(name) == minlen:
                nextset.append(name)
            elif len(name) < minlen:
                nextset = [name]
                minlen = len(name)
        remaining = nextset  # keep the shortest
    if len(remaining) > 1:
        nextset = [name for name in remaining if name.isalnum()]
        if nextset:
            remaining = nextset  # remove non-alphanumerics if possible
    if len(remaining) > 1:
        min_nums = sum([d.isdigit() for d in remaining[0]])
        nextset = [remaining[0]]
        for name in remaining[1:]:
            numdig = sum([d.isdigit() for d in name])
            if numdig == min_nums:
                nextset.append(name)
            elif numdig < min_nums:
                nextset = [name]
                min_nums = numdig
        remaining = nextset  # keep the fewest digits
    if len(remaining) > 1:
        remaining.sort()  # if all else fails, go alphabetical
    chosen = remaining[0]
    if '/' in chosen:
        sys.stderr.write('WARNING: Gene name %s contains illegal character "/" which was replaced with "_"\n' % chosen)
        chosen = chosen.replace('/', '_')  # avoid sub-keying in h5py! Yes, this has happened!
    return chosen

if opts.verbose:
    logprint('Assigning names to transcript families')

if opts.genenames:
    with open(opts.genenames) as infile:
        gene_name_lookup = {x[0]: x[1] for x in [line.strip().split() for line in infile]}
# gene_name_lookup = pd.read_csv(opts.genenames,sep='\t',header=None,names=['tid','tfam']).set_index('tid')['tfam'].to_dict()
else:
    gene_name_lookup = {}

new_tfams = {}
multi_names = defaultdict(lambda: int(1))
for tfam_val in tfams.values():
    geneset = {gene_name_lookup[tid] for tid in tfam_val[0] if tid in gene_name_lookup}
    if not geneset:
        geneset = set(tfam_val[0])  # if no gene names available, just use the tids themselves
    genename = _choose_name(geneset)
    if genename in new_tfams:
        multi_names[genename] += 1
        genename = '%s_%d' % (genename, multi_names[genename])
    new_tfams[genename] = tfam_val
for (genename, num_appearances) in multi_names.items():
    sys.stderr.write('WARNING: Gene name %s appears %d independent times\n' % (genename, num_appearances))

if opts.verbose:
    logprint('Saving results')

with open(outbedname, 'w') as outbed:
    with open(outtxtname, 'w') as outtxt:
        for tfam, (tids, (chrom, strand), genpos) in new_tfams.items():
            outbed.write(SegmentChain(*positionlist_to_segments(chrom, strand, list(genpos)), ID=tfam).as_bed())
            for tid in tids:
                outtxt.write('%s\t%s\n' % (tid, tfam))

if opts.verbose:
    logprint('Tasks complete')
