from unittest import TestCase
from six import assertRaisesRegex
from six.moves import builtins
from os.path import join
from tempfile import mkdtemp
import shutil

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from .mocking import mockOpen

from dark.process import Executor
from dark.reads import Read, Reads

from py3seq import RecombinationAnalysis, readRecombinants
from py3seq.analysis import _OUTPUT_PREFIX, _RECOMBINANTS_HEADER


class TestAnalysis(TestCase):
    """
    Tests for the C{py3seq.RecombinationAnalysis} class.
    """
    @classmethod
    def setUpClass(cls):
        """
        Make a tiny p-value table in a temporary directory.
        """
        e = Executor()
        tmpDir = mkdtemp()
        tableFile = join(tmpDir, 'table')
        e.execute('3seq -g %s 2' % tableFile)
        cls._tableFile = tableFile
        cls._tmpDir = tmpDir

    @classmethod
    def tearDownClass(cls):
        """
        Remove the p-value table temporary directory.
        """
        shutil.rmtree(cls._tmpDir)

    def setUp(self):
        self.ra = RecombinationAnalysis(TestAnalysis._tableFile)

    def tearDown(self):
        if self.ra.tmpDir:
            self.ra.removeOutput()

    def testCheck(self):
        """
        Test that the check function returns a C{CompletedProcess} instance
        with exit status 0.
        """
        result = self.ra.check()
        self.assertEqual(0, result.returncode)

    def testRun(self):
        reads = Reads([
            Read('id1', 'A' * 200 + 'G' * 200),
            Read('id2', 'A' * 400),
            Read('id3', 'G' * 400),
        ])
        result = self.ra.run(reads)
        self.assertEqual(0, result.returncode)

    def testRecombinantFileWithNoRun(self):
        """
        Test that the recombinantFile method raises a RuntimeError if it
        is called before any analysis is done.
        """
        error = '^No analysis has been run yet$'
        assertRaisesRegex(self, RuntimeError, error, self.ra.recombinantFile)

    def testRecombinantFile(self):
        """
        Test that the recombinantFile method produces the expected string.
        """
        reads = Reads([
            Read('id1', 'A' * 200 + 'G' * 200),
            Read('id2', 'A' * 400),
            Read('id3', 'G' * 400),
        ])
        self.ra.run(reads)
        self.assertEqual(
            join(self.ra.tmpDir, _OUTPUT_PREFIX + '.3s.rec'),
            self.ra.recombinantFile())

    def testRemoveOutputWithNoRun(self):
        """
        Test that the removeOutput method raises a RuntimeError if it
        is called before any analysis is done.
        """
        error = '^No analysis has been run yet$'
        assertRaisesRegex(self, RuntimeError, error, self.ra.removeOutput)

    @patch('shutil.rmtree')
    def testRemoveOutput(self, rmtreeMock):
        """
        Test that the removeOutput method is called with the expected
        temporary directory name.
        """
        reads = Reads([
            Read('id1', 'A' * 200 + 'G' * 200),
            Read('id2', 'A' * 400),
            Read('id3', 'G' * 400),
        ])
        self.ra.run(reads)
        self.ra.removeOutput()
        rmtreeMock.assert_called_once_with(self.ra.tmpDir)


class TestReadRecombinants(TestCase):
    """
    Tests for the readRecombinants function.
    """

    def testUnrecognizedHeader(self):
        """
        If an unrecognized header line is passed, readRecombinants must
        raise a ValueError.
        """
        mockOpener = mockOpen(read_data='bad header\n')
        with patch.object(builtins, 'open', mockOpener):
            error = '^Unrecognized header line: bad header$'
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('filename'))

    def testNoRecombinants(self):
        """
        If the recombinants file has only a header, the readRecombinants
        function must not yield any results.
        """
        mockOpener = mockOpen(read_data=('%s\n' % _RECOMBINANTS_HEADER))
        with patch.object(builtins, 'open', mockOpener):
            self.assertEqual([], list(readRecombinants('filename')))

    def testInsufficientFields(self):
        """
        If a line of the recombinants has less than 13 fields, a ValueError
        must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'only\tthree\tfields'
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r'^not enough values to unpack \(expected 13, got 3\)$'
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('filename'))

    def testMNotAnInt(self):
        """
        If a line of the recombinants has an m value that is not an int, a
        ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 six 0 0 1.0 1 3.0 4.0 4.0 5 ...'.replace(' ', '\t'),
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^invalid literal for int\(\) with base 10: 'six'$"
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('filename'))

    def testNNotAnInt(self):
        """
        If a line of the recombinants has an n value that is not an int, a
        ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 six 0 1.0 1 3.0 4.0 4.0 5 ...'.replace(' ', '\t'),
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^invalid literal for int\(\) with base 10: 'six'$"
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('filename'))

    def testKNotAnInt(self):
        """
        If a line of the recombinants has a k value that is not an int, a
        ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 six 1.0 1 3.0 4.0 4.0 5 ...'.replace(' ', '\t'),
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^invalid literal for int\(\) with base 10: 'six'$"
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('filename'))

    def testPNotAFloat(self):
        """
        If a line of the recombinants has a p value that is not a float, a
        ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 0 six 1 3.0 4.0 4.0 5 ...'.replace(' ', '\t'),
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^could not convert string to float: 'six'$"
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('filename'))

    def testHSNotZeroOrOne(self):
        """
        If a line of the recombinants has an hs value that is not '0' or '1', a
        KeyError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 0 3.1 x 3.0 4.0 4.0 5 ...'.replace(' ', '\t'),
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^'x'$"
            assertRaisesRegex(self, KeyError, error, list,
                              readRecombinants('filename'))

    def testLogPNotAFloat(self):
        """
        If a line of the recombinants has a log(p) value that is not a float, a
        ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 0 6 1 six 4.0 4.0 5 ...'.replace(' ', '\t'),
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^could not convert string to float: 'six'$"
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('filename'))

    def testDsPNotAFloat(self):
        """
        If a line of the recombinants has a DS(p) value that is not a float, a
        ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 0 6 1 4.0 4.0 six 5 ...'.replace(' ', '\t'),
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^could not convert string to float: 'six'$"
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('filename'))

    def testMinRecLengthNotAnInt(self):
        """
        If a line of the recombinants has a minimum recombinant length value
        that is not an int, a ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 6 1.0 1 3.0 4.0 4.0 six ...'.replace(' ', '\t'),
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^invalid literal for int\(\) with base 10: 'six'$"
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('filename'))

    def testNoBreakpointsOnLine(self):
        """
        If a line of the recombinants has no breakpoint information, a
        ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 6 1.0 1 3.0 4.0 4.0 6 '.replace(' ', '\t'),
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^No breakpoints found on line 2 of file\.rec$"
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('file.rec'))

    def testBreakpointLeft1NotInteger(self):
        """
        If a line of the recombinants has a first breakpoint index a
        non-integer, a ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 6 1.0 1 3.0 4.0 4.0 6 '.replace(' ', '\t') +
            'a-2 & 3-4',
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^invalid literal for int\(\) with base 10: 'a'$"
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('file.rec'))

    def testBreakpointLeft2NotInteger(self):
        """
        If a line of the recombinants has a second breakpoint index a
        non-integer, a ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 6 1.0 1 3.0 4.0 4.0 6 '.replace(' ', '\t') +
            '2-a & 3-4',
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^invalid literal for int\(\) with base 10: 'a'$"
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('file.rec'))

    def testBreakpointRight1NotInteger(self):
        """
        If a line of the recombinants has a third breakpoint index a
        non-integer, a ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 6 1.0 1 3.0 4.0 4.0 6 '.replace(' ', '\t') +
            '2-3 & a-4',
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^invalid literal for int\(\) with base 10: 'a'$"
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('file.rec'))

    def testBreakpointRight2NotInteger(self):
        """
        If a line of the recombinants has a fourth breakpoint index a
        non-integer, a ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 6 1.0 1 3.0 4.0 4.0 6 '.replace(' ', '\t') +
            '2-3 & 4-a',
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = r"^invalid literal for int\(\) with base 10: 'a'$"
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('file.rec'))

    def testBreakpointIndicesDescendingInPair(self):
        """
        If a line of the recombinants has a pair of breakpoint indices that
        are descending, a ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 6 1.0 1 3.0 4.0 4.0 6 '.replace(' ', '\t') +
            '2-1 & 4-6',
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = (r"^Breakpoints \(2-1 & 4-6\) on line 2 of file.rec "
                     r"do not have non-descending indices$")
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('file.rec'))

    def testBreakpointIndicesDescendingAcrossPair(self):
        """
        If a line of the recombinants has breakpoint indices that
        are descending across a pair of indices, a ValueError must be raised.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 0 6 1.0 1 3.0 4.0 4.0 6 '.replace(' ', '\t') +
            '1-5 & 4-6',
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            error = (r"^Breakpoints \(1-5 & 4-6\) on line 2 of file.rec "
                     r"do not have non-descending indices$")
            assertRaisesRegex(self, ValueError, error, list,
                              readRecombinants('file.rec'))

    def testExpectedRecombinant(self):
        """
        If a line of the recombinants is syntactically correct, the expected
        Recombinant instance must be returned.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 1 6 1.0 1 3.0 5.0 4.0 6 '.replace(' ', '\t') +
            ' 1-3 &  4-6\t10-12 & 50-62',
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            (recombinant,) = list(readRecombinants('file.rec'))

        self.assertEqual('id1', recombinant.pId)
        self.assertEqual('id2', recombinant.qId)
        self.assertEqual('id3', recombinant.recombinantId)
        self.assertEqual(0, recombinant.m)
        self.assertEqual(1, recombinant.n)
        self.assertEqual(6, recombinant.k)
        self.assertEqual(1.0, recombinant.p)
        self.assertEqual(True, recombinant.hs)
        self.assertEqual(3.0, recombinant.logp)
        self.assertEqual(4.0, recombinant.dsP)
        self.assertEqual(6, recombinant.minRecLength)
        self.assertEqual(
            (
                ((1, 3), (4, 6)),
                ((10, 12), (50, 62))
            ),
            recombinant.breakpoints)

    def testTwoExpectedRecombinants(self):
        """
        If two lines of recombinants are syntactically correct, the expected
        Recombinant instances must be returned.
        """
        mockOpener = mockOpen(read_data='\n'.join((
            _RECOMBINANTS_HEADER,
            'id1 id2 id3 0 1 6 1.0 1 3.0 5.0 4.0 6 '.replace(' ', '\t') +
            ' 1-3 &  4-6\t10-12 & 50-62',
            'id4 id5 id6 1 2 7 2.0 0 4.0 6.0 5.0 7 '.replace(' ', '\t') +
            ' 2-4 &  5-7\t11-13 & 51-63',
        )) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            (recombinant1, recombinant2) = list(readRecombinants('file.rec'))

        self.assertEqual('id1', recombinant1.pId)
        self.assertEqual('id2', recombinant1.qId)
        self.assertEqual('id3', recombinant1.recombinantId)
        self.assertEqual(0, recombinant1.m)
        self.assertEqual(1, recombinant1.n)
        self.assertEqual(6, recombinant1.k)
        self.assertEqual(1.0, recombinant1.p)
        self.assertEqual(True, recombinant1.hs)
        self.assertEqual(3.0, recombinant1.logp)
        self.assertEqual(4.0, recombinant1.dsP)
        self.assertEqual(6, recombinant1.minRecLength)
        self.assertEqual(
            (
                ((1, 3), (4, 6)),
                ((10, 12), (50, 62))
            ),
            recombinant1.breakpoints)

        self.assertEqual('id4', recombinant2.pId)
        self.assertEqual('id5', recombinant2.qId)
        self.assertEqual('id6', recombinant2.recombinantId)
        self.assertEqual(1, recombinant2.m)
        self.assertEqual(2, recombinant2.n)
        self.assertEqual(7, recombinant2.k)
        self.assertEqual(2.0, recombinant2.p)
        self.assertEqual(False, recombinant2.hs)
        self.assertEqual(4.0, recombinant2.logp)
        self.assertEqual(5.0, recombinant2.dsP)
        self.assertEqual(7, recombinant2.minRecLength)
        self.assertEqual(
            (
                ((2, 4), (5, 7)),
                ((11, 13), (51, 63))
            ),
            recombinant2.breakpoints)
