## py3seq - a Python wrapper for 3seq recombination detection


Runs under Python 2.7, 3.5, 3.6, and 3.7. [Change log](CHANGELOG.md)
[![Build Status](https://travis-ci.org/acorg/py3seq.svg?branch=master)](https://travis-ci.org/acorg/py3seq)


Here's a Python class that can be used to call the `3seq` recombination
detection algorithm described in
[An exact nonparametric method for inferring mosaic structure in sequence triplets](http://www.genetics.org/content/176/2/1035)
by [Maciej F. Boni](https://www.huck.psu.edu/people/maciej-f-boni),
[David Posada](http://darwin.uvigo.es/dposada/), and
[Marcus W. Feldman](http://www-evo.stanford.edu/marc.html) in
*Genetics*. 2007 Jun; 176(2): 1035–1047. doi:
`10.1534/genetics.106.068874`.

An improvement to the algorithm was described in
[Improved Algorithmic Complexity for the 3SEQ Recombination Detection Algorithm](https://academic.oup.com/mbe/article/35/1/247/4318635)
by
[Ha Minh Lam](https://www.hic-vac.org/members/members-profiles/ha-minh-lam),
[Oliver Ratmann](https://www.imperial.ac.uk/people/oliver.ratmann05), and
[Maciej F. Boni](https://www.huck.psu.edu/people/maciej-f-boni) in
*Molecular Biology and Evolution*, Volume 35, Issue 1, 1 January 2018,
Pages 247–251, doi `10.1093/molbev/msx263`.

## 3seq pre-requisite

The source code for `3seq` can be obtained from
[http://mol.ax/software/3seq](http://mol.ax/software/3seq/). The manual
[is here](http://mol.ax/content/media/2018/02/3seq_manual.20180209.pdf).

You will need to compile the software and then move the resulting `3seq`
exectuable into a directory that's in your shell's `PATH` variable.

### p-value lookup table

You will also need a p-value lookup table for `3seq` to use. You can either
generate one yourself (see steps 3a, 3b, and 4 in
[the manual](http://mol.ax/content/media/2018/02/3seq_manual.20180209.pdf)). Or
else download and uncompress `PVT.3SEQ.2017.700` from
[http://mol.ax/3seq](http://mol.ax/3seq).

## Installing py3seq

```sh
$ pip install py3seq
```

## Using py3seq

```python
from __future__ import print_function

from py3seq import RecombinationAnalysis, readRecombinants

# The p-value lookup table file PVT.3SEQ.2017.700 is available as described above.
analysis = RecombinationAnalysis('PVT.3SEQ.2017.700')

# Run the recombination on a FASTA file (Phylip is also supported, or you
# can pass a Reads instance from the dark-matter package).
analysis.run('filename.fasta')

# The 3seq output files can now be accessed in analysis.tmpDir in case you
# need them. See section 8 of the 3seq manual for their names.

# Process all predicted recombinants.
for recombinant in readRecombinants(analysis.recombinantFile()):
    print('%s is predicted to be a recombinant of %s and %s' %
          (recombinant.recombinantId, recombinant.pId, recombinant.qId))
    # See py3seq for all attributes of the Recombinant class.

# Remove 3seq output files.
analysis.removeOutput()

# See what shell commands were run in performing the analysis and how long they
# took, etc.
print('\n'.join(analysis.executor.log))
```

## Development

```sh
$ git clone https://github/acorg/py3seq
$ cd py3seq
$ pip install -e '.[dev]'   # To install development dependencies.
$ make check
```
