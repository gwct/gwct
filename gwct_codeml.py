#!/usr/bin/python
#############################################################################
#Runs codeml on a single .fa file or a directory full of .fa files.
#
#Dependencies: PAML, newickutils (if you want to prune your tree)
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os, argparse, lib.gwctcore as gc, lib.gwctree as gt
from random import randint

aas = ["G","A","V","L","I","P","F","Y","W","S","T","C","M","N","Q","K","R","H","D","E","X","-"];
nts = ["A","T","C","G","N","-","X"];

############################################
#Function Definitions
############################################
def optParse():
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Runs codeml on a directory full of .fa files. Dependencies: PAML, newickutils (if you want to prune your tree with --prune)");

	parser.add_argument("-i", dest="input", help="Input. A directory containing many FASTA (.fa) files.");
	parser.add_argument("-p", dest="paml_path", help="You must specify the full path to your PAML DIRECTORY here.");
	parser.add_argument("-t", dest="tree_file", help="A user specified tree for codeml to use. If not specified, codeml will infer the tree.", default="");
	parser.add_argument("--prune", dest="prune_opt", help="If not all species present in the tree will be present in each alignment, set this flag to prune the tree for each file.", action="store_true");
	# parser.add_argument("-seqtype", dest="paml_seqtype", help="Enter either 'codon' or 'aa'. Default value is 'codon'.", default='codon');
	parser.add_argument("-v", dest="verbosity", help="An option to control the output printed to the screen. 1: print all codeml output, 0: print only a progress bar. Default: 1", type=int, default=1);
	parser.add_argument("-o", dest="output", help="Desired output directory. If none is entered, will be determined automatically.", default=False);

	args = parser.parse_args();

	if args.input == None or args.paml_path == None:
		sys.exit(gc.errorOut(1, "Both -i and -c must be set."));
	if not os.path.isdir(args.input) or not os.path.isdir(args.paml_path):
		sys.exit(gc.errorOut(2, "Both -i and -c must be valid directory paths!"));
	else:
		args.input = os.path.abspath(args.input);
		args.paml_path = os.path.abspath(args.paml_path);
	if not os.path.isfile(args.tree_file):
		sys.exit(gc.errorOut(3, "-t must be a valid file name."));
	#else:
	#	args.tree_file = os.path.abspath(args.tree_file);
	try:
		td, tree, r = gt.treeParse(open(args.tree_file, "r").read().replace("\n",""));
	except:
		sys.exit(gc.errorOut(4, "-t does not contain a valid Newick string!"));
	# if args.paml_seqtype not in ['codon','aa']:
	# 	sys.exit(gc.errorOut(5, "-seqtype must be either 'codon' or 'aa'."));
	if args.verbosity not in [0,1]:
		sys.exit(gc.errorOut(6, "-v must take values of either 1 or 0"));

	return args.input, args.paml_path, args.tree_file, args.prune_opt, args.verbosity, args.output;

############################################
#Main Block
############################################

indir, path, treefile, prune, v, output = optParse();

starttime = gc.getLogTime();

file_flag = 0;
outdir = gc.defaultOut(output, indir, "-gwct-codeml");
codemldir = os.path.join(outdir, "codeml-out");
ancdir = os.path.join(outdir, "anc-seqs-fa");
logfilename = os.path.join(outdir, "gwct-codeml.log");

pad = 40;
print gc.spacedOut(gc.getTime() + " | Creating main output directory:", pad) + outdir;
os.system("mkdir " + outdir);
gc.printWrite(logfilename, "=======================================================================");
gc.printWrite(logfilename, "\tRunning codeml to reconstruct ancestral sequences");
gc.printWrite(logfilename, "\t\t\t" + gc.getDateTime());
gc.printWrite(logfilename, gc.spacedOut("INPUT    | Running codeml on all files in:", pad) + indir);
gc.printWrite(logfilename, gc.spacedOut("INFO     | PAML path set to:", pad) + path);
if treefile != "":
	gc.printWrite(logfilename, gc.spacedOut("INFO     | Using tree from file:", pad) + treefile);
else:
	gc.printWrite(logfilename, "INFO     | No tree file specified. codeml will infer a tree for each gene.");
if prune:
	gc.printWrite(logfilename, "INFO     | Pruning the tree for each gene.");
# if seqtype == 'codon':
#	gc.printWrite(logfilename, gc.spacedOut("INFO     | Seqtype set to:", pad) + "codons");
# if seqtype == 'aa':
#	gc.printWrite(logfilename, gc.spacedOut("INFO     | Seqtype set to:", pad) + "amino acids (aa)");
gc.printWrite(logfilename, gc.spacedOut("OUTPUT   | Writing output to:", pad) + outdir);
if v == 1:
	gc.printWrite(logfilename, "INFO     | Printing all codeml output to the screen.");
else:
	gc.printWrite(logfilename, "INFO     | Silent mode. Not printing codeml output to the screen.");
gc.printWrite(logfilename, "-------------------------------------");
filelist = os.listdir(indir);
print "** Creating codeml output directory:\t" + codemldir;
os.system("mkdir " + codemldir);
print "** Creating directory to pass ancestral sequences and trees:\t" + ancdir;
os.system("mkdir " + ancdir);

if prune:
	print "** Retrieving tree info...";
	td, tree, r = gt.treeParse(open(treefile, "r").read().replace("\n",""),0);
	tips = [node for node in td if td[node][2] == 'tip'];

gc.printWrite(logfilename, gc.getTime() + " | Starting codeml runs...\n");
if v == 0:
	codeml_logfile = os.path.join(outdir, "codeml.stdout");
ctlfilename = "codeml.ctl";

i, numbars, donepercent, numfiles = 0, 0, [], len(filelist);
fa_skip = [];
for cur_file in filelist:
	if v == 0:
		numbars, donepercent = gc.loadingBar(i, numfiles, donepercent, numbars);
	i += 1;

	if not cur_file.endswith(".fa"):
		fa_skip.append(cur_file);
		continue;
	infilename = os.path.join(indir, cur_file);
	gid = cur_file[:cur_file.index(".")];
	run_outdir = os.path.join(codemldir, gid);

	if not os.path.exists(run_outdir):
		os.system("mkdir " + run_outdir);

	if prune == 1:
		seqs = gc.fastaGetDict(infilename);
		to_prune = [tip for tip in tips if (">" + tip) not in seqs];

		nw_cmd = "nw_prune " + treefile + " ";
		for tip in to_prune:
			nw_cmd += tip + " ";
		nw_cmd += " > pruned.tre";
		os.system(nw_cmd);

	with open(ctlfilename, "w") as ctlfile:
		inline = "seqfile = " + infilename + "\n";
		ctlfile.write(inline);
		if treefile != "":
			if prune:
				treeline = "treefile = pruned.tre\n";
			else:
				treeline = "treefile = " + treefile + "\n";
			ctlfile.write(treeline);
		outline = "outfile = " + os.path.join(run_outdir, gid + ".out") + "\n\n";
		ctlfile.write(outline);

		ctlfile.write("noisy = 3\n");
		ctlfile.write("verbose = 0\n");
		if treefile != "":
			ctlfile.write("runmode = 0\n\n");
		else:
			ctlfile.write("runmode = 2\n\n");

		# if seqtype == 'codon':
		#	ctlfile.write("seqtype = 1\n");
		# elif seqtype == 'aa':
		ctlfile.write("seqtype = 2\n");

		ctlfile.write("CodonFreq = 2\n");
		ctlfile.write("clock = 0\n");
		ctlfile.write("aaDist = 0\n");
		ctlfile.write("aaRatefile = " + os.path.join(path, "dat", "wag.dat") + "\n");
		ctlfile.write("model = 2\n\n");
		ctlfile.write("NSsites = 0\n\n");

		ctlfile.write("icode = 0\n");
		ctlfile.write("fix_kappa = 0\n");
		ctlfile.write("kappa = 3\n");
		ctlfile.write("fix_omega = 0\n");
		ctlfile.write("omega = 1\n\n");

		ctlfile.write("fix_alpha = 1\n");
		ctlfile.write("alpha = 0\n");
		ctlfile.write("Malpha = 0\n");
		ctlfile.write("ncatG = 10\n\n");

		ctlfile.write("getSE = 0\n");
		ctlfile.write("RateAncestor = 1\n");
		ctlfile.write("Small_Diff = .5e-6\n");


	codeml_cmd = os.path.join(path, "bin", "codeml " + ctlfilename);
	if v == 0:
		codeml_cmd += " >> " + codeml_logfile;
		gc.printWrite(logfilename, gc.spacedOut(gc.getTime() + " | codeml Call for " + cur_File + ":", 50) + codeml_cmd, mode=2);
	else:
		gc.printWrite(logfilename, gc.spacedOut(gc.getTime() + " | codeml Call for " + cur_file + ":", 50) + codeml_cmd);
	os.system(codeml_cmd);

	anc_seqfile = os.path.join(ancdir, gid + "-anc.fa");
	anc_probfile = os.path.join(ancdir, gid + "-ancprobs.fa");
	anc_treefile = os.path.join(ancdir, gid + "-anc.tre");
	rstfile = open("rst", "r");
	rstlines = rstfile.readlines();
	rstfile.close();

	node_list = [];
	for k in xrange(len(rstlines)):
		if rstlines[k] == "tree with node labels for Rod Page's TreeView\n":
			anctree = rstlines[k+1].replace(" ","");
			atfile = open(anc_treefile, "w");
			atfile.write(anctree);
			atfile.close();

		if rstlines[k] == "List of extant and reconstructed sequences\n":
			asfile = open(anc_seqfile, "w");

			j = 4;
			while rstlines[k+j].find("Overall accuracy of the") == -1:
				if rstlines[k+j] != "\n":
					tmpline = rstlines[k+j].replace("\n", "");
					tmpline = tmpline.split("  ");
					tmpline = filter(None, tmpline);
					node_list.append(tmpline[0]);

					title = ">" + tmpline[0] + "\n";
					asfile.write(title);
					seq = tmpline[1].replace(" ","") + "\n";
					asfile.write(seq);
				j = j + 1;
			asfile.close();
			break;

	curseqs = {};
	for n in node_list:
		curseqs[n] = "";

	for k in xrange(len(rstlines)):
		if rstlines[k] == "Prob of best state at each node, listed by site\n":
			j = 4;
			while rstlines[k+j] != "\n":
				tmpline = filter(None, rstlines[k+j].replace("\n","").split("  "));

				cur_states = tmpline[2].split(": ");
				extant = list(cur_states[0].replace(" ",""));
				ancs = filter(None, cur_states[1].split(" "));
				all_states = extant + ancs;

				for n in xrange(len(node_list)):
					cur_spec = node_list[n];
					cur_aa = all_states[n];
					if cur_aa.find("(") != -1 and cur_aa.find(")") != -1:
						cur_aa = cur_aa.replace("(","_").replace(")","") + " ";
					curseqs[cur_spec] = curseqs[cur_spec] + cur_aa;
				j = j + 1;

	asfile = open(anc_probfile, "w");
	for seq in curseqs:
		title = ">" + seq + "\n";
		asfile.write(title);
		for base in curseqs[seq]:
			asfile.write(base);
		asfile.write("\n");
	asfile.close();


	newfileList = os.listdir(os.getcwd());
	for neweach in newfileList:
		if neweach in ["2NG.dN","2NG.dS","2NG.t","codeml.ctl","lnf","rst","rst1","rub","pruned.tre"]:
			mv_cmd = "mv " + neweach + " " + run_outdir;
			os.system(mv_cmd);


if v == 0 and fileflag == 0:
	pstring = "100.0% complete.\n";
	sys.stderr.write('\b' * len(pstring) + pstring);
gc.printWrite(logfilename, gc.getTime() + " | Done!");
gc.printWrite(logfilename, gc.getLogTime());
gc.printWrite(logfilename, "=======================================================================");
