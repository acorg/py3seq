# Note that the version string must have the following format, otherwise it
# will not be found by the version() function in ../setup.py
__version__ = '1.1.1'

from .analysis import RecombinationAnalysis, readRecombinants

# Keep Python linters quiet.
_ = (RecombinationAnalysis, readRecombinants)
