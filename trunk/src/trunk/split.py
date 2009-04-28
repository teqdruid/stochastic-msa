#!/usr/bin/python

import sys

files = dict()
def getFile (seg):
    if not files.has_key(seg):
        files[seg] = file(sys.argv[1] + "seg-" + seg + ".fasta", "w")
    return files[seg]

cfile = None
line = sys.stdin.readline()
while line != "":
    if line[0] == '>':
        seg = "other"
        if (line.find("segment") != -1):
            seg = line.split("segment")[1][:2].lstrip()
        else:
            print line
        cfile = getFile(seg)
    
    if cfile != None:
        cfile.write(line)

    line = sys.stdin.readline()
        

for v in files.values():
    v.close()
