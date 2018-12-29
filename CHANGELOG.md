## 1.1.4 2018-12-29

Allow passing a string value for `t` to the `run` method.

## 1.1.3 2018-12-13

Added a -t argument to `RecombinationAnalysis.run()`.

## 1.1.2 2018-11-23

Added note about license on underlying `3seq` code to README.

## 1.1.1 2018-11-18

Fixed variable name typo in `README.md` example Python.

## 1.1.0 2018-11-18

Use the `3seq` option `-ptable` when running with `-full` to set the
p-value table that should be used. So passing a p-value table file to
`RecombinationAnalysis` (now mandatory) now works correctly. *Note* though
that this will result in the file name being recorded in your `3seq` config
file, replacing whatever value you may have had there previously while
doing other interactive work with `3seq`.

## 1.0.3 2018-11-18

TravisCI tweaks. Changed default location of p-value lookup table
file. This is only used in running tests at the moment.  Changes will be
needed once `3seq` can be given a p-value table file on the command line.

## 1.0.2 2018-11-18

Added TravisCI support.

## 1.0.1 2018-11-18

Initial release.
