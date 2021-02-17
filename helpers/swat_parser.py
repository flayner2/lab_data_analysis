"""
A Python script to parse the results of an alignment produced by running `swat`.
It searches the given folder for a .alignment and a .allscores file. A .alignment
file can be produced by catching `swat`'s printed results to the STDOUT to a file.
A score value may be provided if you also want to use this to find the starting and
ending positions of an alignment.
"""
