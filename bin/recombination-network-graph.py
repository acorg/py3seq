#!/usr/bin/env python

from __future__ import print_function, division

import sys
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt

from dark.dna import _pp
from dark.fasta import FastaReads

from py3seq import (
    RecombinationAnalysis, readRecombinants, triplet, informativeSites)


def identity(dnaMatch, read1, read2, indent=''):
    """
    Return the percent identity between two reads as a string.

    @param dnaMatch: A C{dict} returned by C{compareDNAReads}.
    @param read1: A C{Read} instance or an instance of one of its subclasses.
    @param read2: A C{Read} instance or an instance of one of its subclasses.
    @param indent: A C{str} to indent all returned lines with.
    @return: A C{str} describing the match.
    """
    return _pp('%sExact matches' % indent,
               dnaMatch['match']['identicalMatchCount'],
               len(read1), len(read2))


def tooClose(offsets, seen, fraction):
    """
    Detect whether a breakpoint pair is sufficiently different from *all*
    members of a list of previously analyzed breakpoints. If either of the
    offsets is sufficiently different from all previously seen breakpoints,
    return C{True}.

    @param offsets: A 2-tuple of C{int} start and end offsets.
    @param seen: An iterable of 2-tuples of C{int} start and end offsets.
    @param fraction: The C{float} fraction that is considered too close (i.e.,
        not sufficiently different).
    @return: The 2-tuple that disqualifies the passed offsets, or C{None}.
    """
    for seenStart, seenEnd in seen:
        startMargin = fraction * seenStart
        endMargin = fraction * seenEnd
        # If the passed offsets are both too close to the offets of this
        # previously seen pair, the passed offsets are not sufficiently far
        # away from all others.
        if (abs(offsets[0] - seenStart) <= startMargin and
                abs(offsets[1] - seenEnd) <= endMargin):
            return (seenStart, seenEnd)


parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description=('Print simple recombination analysis results '
                 'for an alignment file'))

parser.add_argument(
    '--alignmentFile', required=True,
    help='The FASTA or Phylip alignment file.')

parser.add_argument(
    '--pValueFile', required=True,
    help='The 3seq p-value lookup table file.')

parser.add_argument(
    '--verbose', default=False, action='store_true',
    help='Print information about the ongoing analysis to standard error')

parser.add_argument(
    '--keep', default=False, action='store_true',
    help='Do not remove the 3seq recombinant file. Instead print its path.')

parser.add_argument(
    '--show', default=False, action='store_true',
    help='Show the network graph')

parser.add_argument(
    '--imageFile',
    help='The file name to save the recombination graph to.')

parser.add_argument(
    '--tValue', '-t', default='0.05',
    help='The -t value to pass to 3seq')

parser.add_argument(
    '--strict', default=False, action='store_true',
    help='If given, do not allow ambiguous nucleotide symbols to match')

parser.add_argument(
    '--breakpointSimilarityThreshold', type=float, default=0.95,
    help=('If both offsets of a breakpoint pair are not sufficiently '
          'different from previously printed breakpoints, the breakpoint. '
          'The breakpoint will not be printed. This [0..1] argument controls '
          'how numerically close breakpoints will be to be considered too '
          'close.'))

args = parser.parse_args()

reads = list(FastaReads(args.alignmentFile))

if not reads:
    print('No reads found in %s!' % args.alignmentFile)
    sys.exit(1)

readsDict = dict((read.id, read) for read in reads)
sequenceLength = len(reads[0])
matchAmbiguous = not args.strict
analysis = RecombinationAnalysis(args.pValueFile)
if args.verbose:
    print('Starting 3seq analysis...', file=sys.stderr)
analysis.run(args.alignmentFile, t=args.tValue)
if args.verbose:
    print('Done.', file=sys.stderr)
if args.show or args.imageFile:
    import networkx as nx
    from networkx.drawing.layout import circular_layout
    G = nx.DiGraph(directed=True)
seen = set()
parentPairs = defaultdict(list)
count = 0

for count, recombinant in enumerate(
        readRecombinants(analysis.recombinantFile()), start=1):
    if count > 1:
        print()
    print('%s may be a recombinant of %s and %s.' %
          (recombinant.recombinantId, recombinant.pId, recombinant.qId))
    p = readsDict[recombinant.pId]
    q = readsDict[recombinant.qId]
    child = readsDict[recombinant.recombinantId]
    iSites = informativeSites(p, q, child)
    print('  There were %d informative sites overall.' % len(iSites))
    parents = tuple(sorted([recombinant.pId, recombinant.qId]))
    if args.show or args.imageFile:
        G.add_edges_from(((recombinant.pId, recombinant.recombinantId),
                          (recombinant.qId, recombinant.recombinantId)))
    key = (parents, recombinant.recombinantId)
    if key in seen:
        print('  ALREADY SEEN!')
    else:
        seen.add(key)

    parentPairs[parents].append(recombinant)

    offsetsPrinted = []

    for breakpointCount, breakpoint in enumerate(
            recombinant.breakpoints, start=1):

        start = (breakpoint[0][0] + breakpoint[0][1]) >> 1
        end = (breakpoint[1][0] + breakpoint[1][1]) >> 1

        print('  Breakpoint %d: %d-%d & %d-%d (using start %d end %d)' % (
            breakpointCount,
            breakpoint[0][0], breakpoint[0][1],
            breakpoint[1][0], breakpoint[1][1],
            start, end))

        closest = tooClose((start, end), offsetsPrinted,
                           args.breakpointSimilarityThreshold)
        if closest:
            print('    Not printing (too close to previous breakpoint '
                  '(%d, %d)).' % closest)
            continue

        offsetsPrinted.append((start, end))

        breakpoints = []
        if start:
            breakpoints.append(start)

        if end == sequenceLength:
            breakpoints.append(end)
        else:
            breakpoints.extend((end, sequenceLength))

        a = triplet(p, q, child, breakpoints=breakpoints,
                    matchAmbiguous=matchAmbiguous)

        for breakpointRegionCount, ba in enumerate(a['breakpointAnalysis'],
                                                   start=1):
            print('    Region %d: start = %d, end = %d, '
                  'length = %d' %
                  (breakpointRegionCount, ba['start'], ba['end'],
                   ba['length']))
            print('      %d informative sites in region (m = %d, n = %d).' %
                  (len(ba['informativeSites']), ba['m'], ba['n']))

            # The informative sites in the region.
            print('      Informative sites: p -> c:', identity(
                ba['informative']['p child'],
                ba['pRegionInformative'], ba['childRegionInformative']))
            print('                         q -> c:', identity(
                ba['informative']['q child'],
                ba['qRegionInformative'], ba['childRegionInformative']))

            # All sites in the region.
            print('      All sites:         p -> q:', identity(
                ba['all']['p q'], ba['pRegion'], ba['qRegion']))
            print('                         p -> c:', identity(
                ba['all']['p child'], ba['pRegion'], ba['childRegion']))
            print('                         q -> c:', identity(
                ba['all']['q child'], ba['qRegion'], ba['childRegion']))

if args.show or args.imageFile:
    nx.draw_networkx(G, pos=circular_layout(G), node_color='blue', width=1,
                     node_size=100, arrowstyle='-|>', arrowsize=12,
                     font_size=8)
    plt.axis('off')
    plt.title('%d recombinant%s (t = %s)' % (
        count, '' if count == 1 else 's', args.tValue))
    if args.show:
        plt.show()
    if args.imageFile:
        plt.savefig(args.imageFile)

if args.keep:
    print('Analysis recombinant file saved to', analysis.recombinantFile())
else:
    analysis.removeOutput()

print('%s potential recombinant%s found.' % (count, '' if count == 1 else 's'))
