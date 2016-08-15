import sys, os

s = 11;
l = 55;

resultsdir = "/Users/Gregg/Desktop/marine_2/seq_files/june/uniq_reps_proper/";

outfilename = resultsdir + "mm_reps_proper_" + str(s) + "-" + str(l-1) + "_counts.txt";
outfile = open(outfilename, "w");
outfile.write("");
counts = {};

for x in xrange(s,l):
	print x;
	counts[x] = [];
	indir = resultsdir + str(x) + "spec_100reps/";
	#indir = resultsdir + str(x) + "spec/";
	filelist = os.listdir(indir);
	for each in filelist:
		infilename = indir + each;
		for line in open(infilename):
			if line.find("# Total sites found: ") != -1:
				line = line.replace("\n","");
				count = line[line.index(":")+2:];
			if line.find("# Total genes checked: ") != -1:
				line = line.replace("\n","");
				total = line[line.index(":")+2:];
		percent = float(count)/float(total);
#		counts[x].append(percent);

		counts[x].append(count);


for x in xrange(s,l):
	outfile.write(str(x) + " ");
outfile.write("\n");

for y in xrange(0,100):
	for x in xrange(s,l):
		#print x;
		outline = counts[x][y], " ";
		outfile.write(str(counts[x][y]) + " ");
	outfile.write("\n");
outfile.close();
