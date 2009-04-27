#!/usr/bin/python

maxLineLen = 80

import sys
from cStringIO import StringIO

# First, read in sequences

outfile = sys.argv[1]

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

out = file(outfile, "w")

print >> out, "!!NA_MULTIPLE_ALIGNMENT 1.0"
print >> out
print >> out, outfile, ".."
print >> out
for (k, v) in seqs.items():
    if (k != ""):
        s = v.getvalue()
        print >> out, ("Name: " + nameFormat + " Len: %d Check: %d Weight: 1.0") % (k, len(s), 0)
    else:
        seqs.pop(k)

print >> out
print >> out, "//"
print >> out

while has_more():
    for (k, v) in seqs.items():
        if (k != ""):
            out.write(nameFormat % (k))
            if v != None:
                chars = v.read(charsPerLine)
                out.write(chars)
                if len(chars) < charsPerLine:
                    seqs[k] = None
            out.write("\n")
    out.write("\n")

out.close()
