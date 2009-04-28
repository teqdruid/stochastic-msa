#!/usr/bin/python

import sys

if len(sys.argv) != 2:
    print "Usage: %s [test aln]" % sys.argv[0]
    sys.exit(1)

tseqs = set()

ta = file(sys.argv[1], "r")

line = ta.readline()
while line.find("//") == -1:
    if line.find("Name:") != -1:
        name = line.split("Name:")[1].lstrip().split(" ")[0]
        tseqs.add(name)
    line = ta.readline()
    

remove = set()

line = sys.stdin.readline()
while line.find("//") == -1:
    if line.find("Name:") != -1:
        name = line.split("Name:")[1].lstrip().split(" ")[0]
        if not name in tseqs:
            remove.add(name)
        else:
           sys.stdout.write(line)
    else:
        sys.stdout.write(line)
    line = sys.stdin.readline()


for v in remove:
    print >> sys.stderr, "Removing ", v

def skip_line(l):
    for v in remove:
        if line.startswith(v):
            return True
    return False

while line != "":
    if not skip_line(line):
        sys.stdout.write(line)
    line = sys.stdin.readline()
