#!/usr/bin/python
#############################################################################
# Script to count unique substitutions among a set of species.
#
# Gregg Thomas
# June 2016
#############################################################################

import sys, os, argparse, lib.gwctcore as gwctcore

############################################
#Function Definitions
############################################
def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Unique substitution counting.");

	parser.add_argument("-i", dest="input", help="Input directory containing a single run_codeml output directory or codeml_combined directory.");
	parser.add_argument("-t", dest="target_specs", help="A list contained in quotes. Spaces separate groups, commas separate species within groups. E.g: \"triMan1 lepWed1,odoRosDiv1 orcOrc1,turTru2\"");
	parser.add_argument("-o", dest="output_suffix", help="The SUFFIX of the directory name to be created by the script. Something descriptive for your run.", default="");
	
	args = parser.parse_args();

	if errorflag == 0:

		if args.input == None:
			gwctcore.errorOut(1, "-i must be defined");
			optParse(1);

		if args.target_specs == None:
			gwctcore.errorOut(2, "-t must be defined");
			optParse(1);

		return args.input, args.target_specs, args.output_suffix;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

############################################
#Main Block
############################################

starttime = gwctcore.getLogTime();
indir, targets, suffix = optParse(0);
targets = targets.split(",");

indir, script_outdir = gwctcore.getOutdir(indir, "gwct", starttime, suffix);

print "# " + gwctcore.getTime() + " | Creating main output directory:\t" + script_outdir;
os.system("mkdir " + script_outdir);

logfilename = os.path.join(script_outdir, "gwct.log");
logfile = open(logfilename, "w");
logfile.write("");
logfile.close();
uniq_header = "#\tChromosome\tGeneID\tAlignLen\tPosition\tBackgroundAlleles\tTargetAlleles\n";
l = 1;
sp = 65;

gwctcore.logCheck(l, logfilename, "# =======================================================================");
gwctcore.logCheck(l, logfilename, "#\tGWCT: Genome-Wide Convergence Tester -- Unique substitutions");
gwctcore.logCheck(l, logfilename, "#\t\t\t" + gwctcore.getDateTime());
gwctcore.logCheck(l, logfilename, gwctcore.spacedOut("# INPUT  | Using sequences in:", sp) + indir);
gwctcore.logCheck(l, logfilename, gwctcore.spacedOut("# INPUT  | Target species:", sp) + ",".join(targets));
gwctcore.logCheck(l, logfilename, gwctcore.spacedOut("# OUTPUT | Output directory created within input directory:", sp) + script_outdir);
gwctcore.logCheck(l, logfilename, "# ---------------------------------------------");
gwctcore.logCheck(2, logfilename, uniq_header);
################################
print "# " + gwctcore.getTime() + " | Counting unique substitutions....";

filelist = os.listdir(indir);
numfiles = len(filelist);
numbars = 0;
donepercent = [];
i = 0;

numsubs = 0;
numgenes = 0;

for each in filelist:
	if each.find(".fa") == -1:
		continue;

	gflag = 0;

	numbars, donepercent = gwctcore.loadingBar(i, numfiles, donepercent, numbars);
	i = i + 1;

	infilename = os.path.join(indir, each);
	inseqs = gwctcore.fastaGetDict(infilename);

	# print infilename;

	seqlen = len(inseqs[">" + targets[0]]);
	fid = each[:each.index(".fa")];
	gid = "_".join(fid.split("_")[:2]);
	chrome = fid[fid.find("chr"):fid.find("chr")+4]

	for x in range(len(inseqs[">" + targets[0]])):
		site = [];
		target_alleles = [];
		bg_alleles = [];

		for title in inseqs:
			site.append(inseqs[title][x]);
			if title[1:] in targets:
				target_alleles.append(inseqs[title][x]);
			else:
				bg_alleles.append(inseqs[title][x]);

		if target_alleles.count(target_alleles[0]) == len(target_alleles) and target_alleles[0] not in bg_alleles and "-" not in site and "X" not in site:
			numsubs += 1;
			if gflag == 0:
				gflag = 1;
				numgenes += 1;
			outline = str(numsubs) + "\t" + chrome + "\t" + gid + "\t" + str(seqlen) + "\t" + str(x+1) + "\t" + "".join(bg_alleles) + "\t" + "".join(target_alleles);
			gwctcore.logCheck(2, logfilename, outline);
			#sys.exit();

pstring = "100.0% complete.";
sys.stderr.write('\b' * len(pstring) + pstring);
print
gwctcore.logCheck(l, logfilename, "# ---------------------------------------------");
print len(filelist), "\t", numsubs, "\t", numgenes
gwctcore.logCheck(l, logfilename, "# " + gwctcore.getLogTime() + " | Done!");
gwctcore.logCheck(l, logfilename, "# =======================================================================");













