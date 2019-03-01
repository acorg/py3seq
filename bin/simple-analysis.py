#!/usr/bin/env python

from __future__ import print_function

import sys
import argparse
from collections import defaultdict

from py3seq import RecombinationAnalysis, readRecombinants


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
    '--tValue', '-t',
    help='The -t value to pass to 3seq')

args = parser.parse_args()

analysis = RecombinationAnalysis(args.pValueFile)

if args.verbose:
    print('Starting 3seq analysis...', file=sys.stderr)

analysis.run(args.alignmentFile, t=args.tValue)

if args.verbose:
    print('Done.', file=sys.stderr)

seen = set()
parentPairs = defaultdict(list)

count = 0

for count, recombinant in enumerate(
        readRecombinants(analysis.recombinantFile()), start=1):
    print('%s may be a recombinant of %s and %s' %
          (recombinant.recombinantId, recombinant.pId, recombinant.qId))
    parents = tuple(sorted([recombinant.pId, recombinant.qId]))
    key = (parents, recombinant.recombinantId)
    if key in seen:
        print('  ALREADY SEEN!')
    else:
        seen.add(key)

    parentPairs[parents].append(recombinant)

    for breakpoint in recombinant.breakpoints:
        print('  %d-%d & %d-%d' % (breakpoint[0][0], breakpoint[0][1],
                                   breakpoint[1][0], breakpoint[1][1]))

analysis.removeOutput()

print('%s potential recombinant%s found.' % (count, '' if count == 1 else 's'))

if count:
    print('Summary of %d parent pair%s:' %
          (len(parentPairs), '' if len(parentPairs) == 1 else 's'))

    for (pId, qId), recombinants in parentPairs.items():
        print('%d possible recombinant%s from parents: %s and %s' %
              (len(recombinants), '' if len(recombinants) == 1 else 's',
               pId, qId))
        for recombinant in recombinants:
            print('  Child %s' % recombinant.recombinantId)
            for breakpoint in recombinant.breakpoints:
                print('    %d-%s & %d-%d' %
                      (breakpoint[0][0], breakpoint[0][1],
                       breakpoint[1][0], breakpoint[1][1]))
