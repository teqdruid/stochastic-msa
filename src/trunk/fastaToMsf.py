#!/usr/bin/python

maxLineLen = 80

import sys
from cStringIO import StringIO

# First, read in sequences

seqs = dict()
currSeq = StringIO()
currSeqName = ""
maxNameLen = 0

def flush():
    if currSeq != "":
        seqs[currSeqName] = StringIO(currSeq.getvalue())

line = sys.stdin.readline()
while line != "":
    if line[0] == ">":
        flush()
        currSeq = StringIO()
        currSeqName = line[1:(line[1:].find(' '))+1]
        if len(currSeqName) > maxNameLen:
            maxNameLen = len(currSeqName)
    else:
        currSeq.write(line.rstrip())

    line = sys.stdin.readline()

flush()

def has_more():
    for (k, v) in seqs.items():
        if v != None:
            return True
    return False

charsPerLine = maxLineLen - (maxNameLen + 2)
while charsPerLine < 10:
    maxLineLen = maxLineLen + 10
    charsPerLine = maxLineLen - (maxNameLen + 2)

nameFormat = "%-" + str(maxNameLen + 2) + "s"

while has_more():
    for (k, v) in seqs.items():
        sys.stdout.write(nameFormat % (k))
        if v != None:
            chars = v.read(charsPerLine)
            sys.stdout.write(chars)
            if len(chars) < charsPerLine:
                seqs[k] = None
        sys.stdout.write("\n")
    sys.stdout.write("\n")
