#!/Users/aholik/tools/miniconda2/bin/python
import sys

scorefile = str(sys.argv[1])

import pyBigWig
bw = pyBigWig.open(scorefile)
scores = open("scores.txt", "w")
for line in open("positions.bed"):
	cols = line.strip().split()
	vals = bw.values(cols[0], int(cols[1])-1, int(cols[2]))
	print >> scores, vals
scores.close()
bw.close()




