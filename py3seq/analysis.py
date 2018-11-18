# xfrom __future__ import print_function

from os.path import join
from tempfile import mkdtemp
import shutil
import six

from dark.process import Executor

import py3seq

_DEFAULT_PVALUE_FILE = join(py3seq.__file__, 'data', 'PVT.3SEQ.2017.700')
_OUTPUT_PREFIX = 'output'

_RECOMBINANTS_HEADER = '\t'.join(
    ('P_ACCNUM Q_ACCNUM C_ACCNUM m n k p HS? log(p) DS(p) DS(p) '
     'min_rec_length breakpoints').split())


class RecombinationAnalysis(object):
    """
    Perform a 3seq recombination analysis.

    @param pValueFile: The C{str} file name containing precomputed p-values
        (as generated by 3seq -g).
    @param dryRun: If C{True} do not execute any 3seq commands, just log what
        would have been run (see self.executor.log for details).
    """

    def __init__(self, pValueFile=_DEFAULT_PVALUE_FILE, dryRun=False):
        self.pValueFile = pValueFile
        self.tmpDir = None
        self.executor = Executor(dryRun=dryRun)

    def check(self):
        """
        Use the -check function to ensure a correct p-value table can be
        checked.

        @return: A C{subprocess.CompletedProcess} instance.
        """
        return self.executor.execute('3seq -check %s' % self.pValueFile)

    def run(self, reads):
        """
        Run 3seq on some reads. Set self.tmpDir as a side-effect.

        @param reads: Either a C{dark.reads.Reads} instance or a C{str}
            filename.
        @return: A C{subprocess.CompletedProcess} instance.
        """
        self.tmpDir = mkdtemp()

        if isinstance(reads, six.string_types):
            inputFile = reads
        else:
            inputFile = join(self.tmpDir, 'input.fasta')
            reads.save(inputFile, format_='fasta')

        return self.executor.execute(
            'echo y | 3seq -f "%s" -id "%s"' %
            (inputFile, join(self.tmpDir, _OUTPUT_PREFIX)))

    def recombinantFile(self):
        """
        Get the name of the main 3seq recombination output file.

        @raise RuntimeError: If no analysis has been run.
        @return: A C{str} path to the output file.
        """
        if self.tmpDir is None:
            raise RuntimeError('No analysis has been run yet')
        else:
            # The string in the following is always used by 3seq.
            return join(self.tmpDir, _OUTPUT_PREFIX + '.3s.rec')

    def removeOutput(self):
        """
        Remove 3seq output files.

        @raise RuntimeError: if no analysis has been run.
        """
        if self.tmpDir is None:
            raise RuntimeError('No analysis has been run yet')
        else:
            shutil.rmtree(self.tmpDir)


class Recombinant(object):
    """
    Hold information about a recombinant found by 3seq. See section 8 of
    http://mol.ax/content/media/2018/02/3seq_manual.20180209.pdf for a fuller
    description of these fields.

    @param pId: The C{str} sequence id of parent p.
    @param qId: The C{str} sequence id of parent q.
    @param recombinantId: The C{str} sequence id of the child (the
        recombinant).
    @param m: The C{int} number of 'up' steps in the random walk view
        of the sequence of matches between the child and the closest parent.
        These are the (parent) p informative sites.
    @param n: The C{int} number of 'down' steps in the random walk view
        of the sequence of matches between the child and the closest parent.
        These are the (parent) q informative sites.
    @param k: The C{int} maximum descent observed for the triplet.
    @param p: The C{float} uncorrected p-value for this triplet.
    @param hs: C{True} if the Hogan-Siegmund approximation was used for this
        triple else C{False}.
    @param logp: The C{float} log base 10 of the p value.
    @param dsP: The C{float} Dunn-Sidak correction of p.
    @param minRecLength: The C{int} minimum length of the recombinant segments.
        I.e. the length of the shorter of the two recombinant segments.
    @param breakpoints: A C{tuple} of C{tuple}s of C{tuple}s of two C{int}s.
        E.g.,
            (
              ((23, 24), (30, 32)),  # Recombination 1 breakpoints.
              ((88, 88), (95, 97)),  # Recombination 2 breakpoints.
            )
        where each second-level C{tuple} is a pair of predicted breakpoints,
        in which each of the pair is a C{tuple} containing two C{int}s giving
        the boundary of the left or right side of the breakpoint. These are
        all potential breakpoints that minimize expression (4) in the Boni et
        al. (2007) Genetics paper (see ../README.md).
    """

    def __init__(self, pId, qId, recombinantId, m, n, k, p, hs, logp, dsP,
                 minRecLength, breakpoints):
        self.pId = pId
        self.qId = qId
        self.recombinantId = recombinantId
        self.m = m
        self.n = n
        self.k = k
        self.p = p
        self.hs = hs
        self.logp = logp
        self.dsP = dsP
        self.minRecLength = minRecLength
        self.breakpoints = breakpoints


def readRecombinants(filename):
    """
    Read a 3seq recombinant file. The output format is described at
    http://mol.ax/content/media/2018/02/3seq_manual.20180209.pdf

    @raise ValueError: If 1) the input file has an unrecognized header, 2) a
        set of breakpoint indices is not non-descending, 3) no breakpoints
        are found on an imput line, or 4) an input line does not have
        sufficient fields.
    @raise KeyError: If C{hs} is not '0' or '1'.
    @return: A generator that yields C{Recombinant} instances.
    """
    with open(filename) as fp:
        header = fp.readline()[:-1]
        if header != _RECOMBINANTS_HEADER:
            raise ValueError('Unrecognized header line: %s' % header)

        hsDict = {'0': False, '1': True}

        for lineNumber, line in enumerate(fp, start=2):
            # The 3s.rec output file has a minimum of 13 columns.
            (pId, qId, cId, m, n, k, p, hs, logp, _, dsP,
             minRecLength, breakpointsStr) = line.split('\t', maxsplit=12)

            # Explicitly convert to the types we need one by one. This will
            # cause a more easily locatable error than if we do them all at
            # once when createing the Recombinant instance below. The dict
            # in the hs conversion is to force a KeyError if hs is not '0'
            # or '1'.
            m = int(m)
            n = int(n)
            k = int(k)
            p = float(p)
            hs = hsDict[hs]
            logp = float(logp)
            dsP = float(dsP)
            minRecLength = int(minRecLength)

            # Extract all breakpoints pairs. These are separated by TAB and
            # have their offsets justified by spaces.
            breakpoints = []
            for breakpoint in breakpointsStr.split('\t'):
                breakpoint = breakpoint.strip()
                if breakpoint:
                    breakpoints.append(breakpoint)

            breakpointTuples = []
            if breakpoints:
                # Breakpoint pairs are split with an ampersand.
                for breakpoint in breakpoints:
                    offsetsLeft, offsetsRight = map(
                        str.strip, breakpoint.split('&'))
                    # And each of these is an integer range split by '-'.
                    left1, left2 = map(int, offsetsLeft.split('-'))
                    right1, right2 = map(int, offsetsRight.split('-'))
                    # Sanity check
                    if left1 <= left2 < right1 <= right2:
                        breakpointTuples.append(
                            ((left1, left2), (right1, right2)))
                    else:
                        raise ValueError(
                            'Breakpoints (%s) on line %d of %s do not have '
                            'non-descending indices' %
                            (breakpoint, lineNumber, filename))
            else:
                raise ValueError('No breakpoints found on line %d of %s' %
                                 (lineNumber, filename))

            yield Recombinant(
                pId, qId, cId, m, n, k, p, hs, logp, dsP, minRecLength,
                tuple(breakpointTuples))
