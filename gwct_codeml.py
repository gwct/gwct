#!/usr/bin/python
#############################################################################
#Runs codeml on a single .fa file or a directory full of .fa files.
#
#Dependencies: PAML, newickutils (if you want to prune your tree)
#
#Gregg Thomas, Summer 2015
#############################################################################

import sys, os, argparse, lib.gwctcore as gwctcore, lib.gwctree as gwctree
from random import randint

aas = ["G","A","V","L","I","P","F","Y","W","S","T","C","M","N","Q","K","R","H","D","E","X","-"];
nts = ["A","T","C","G","N","-","X"];

############################################
#Function Definitions
############################################
def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Runs codeml on a single .fa file or a directory full of .fa files. Dependencies: PAML, newickutils (if you want to prune your tree)");

	parser.add_argument("-i", dest="input", help="Input. Either a directory containing many FASTA files or a single FASTA file.");
	parser.add_argument("-c", dest="paml_path", help="You must specify the full path to your PAML DIRECTORY here.");
	parser.add_argument("-t", dest="tree_file", help="A user specified tree for codeml to use. If not specified, codeml will infer the tree.", default="");
	parser.add_argument("-p", dest="prune_opt", help="If not all species present in the tree will be present in each alignment, set this to 1 to prune the tree for each file. Default: 0", type=int, default=0);
	parser.add_argument("-s", dest="paml_seqtype", help="The seqtype for codeml to use. 1 (default): Codons; 2: Amino Acids", type=int, default=1);
	parser.add_argument("-b", dest="branch_site", help="Specifies the type of run for PAML's branch site test. 0 (default): Do not do branch site test; 1: Do the null model of the branch site test (model=2, NSsite=2, fix_omega=1, omega=1); 2: Do the alternate model of the branch site test (model=2, NSsite=2, fix_omega=0, omega=1). A branch must be specified in your tree file.", type=int, default=0);
	parser.add_argument("-a", dest="anc_opt", help="Option to tell PAML to do ancestral reconstruction (1) or not (0). Default: 0.", type=int, default=0);
	parser.add_argument("-v", dest="verbosity", help="An option to control the output printed to the screen. 1: print all codeml output, 0: print only a progress bar. Default: 1", type=int, default=1);
	parser.add_argument("-l", dest="log_opt", help="A boolean option to tell the script whether to create a logfile (1) or not (0). Default: 1", type=int, default=1);
	parser.add_argument("-x", dest="logdir_suffix", help="A string to add on to the end of the output directory.");

	args = parser.parse_args();

	if errorflag == 0:

		if args.input == None or args.paml_path == None:
			gwctcore.errorOut(1, "Both -i and -c must be set");
			optParse(1);

		if not os.path.isfile(args.tree_file):
			gwctcore.errorOut(2, "-t must be a valid tree file name");
			optParse(1);

		if args.prune_opt not in [0,1]:
			gwctcore.errorOut(3, "-p must take values of either 1 or 0");
			optParse(1);

		if args.prune_opt == 1 and args.tree_file == "":
			gwctcore.errorOut(4, "With -p set to 1 a tree file must be specified");
			optParse(1);

		if args.paml_seqtype not in [1,2]:
			gwctcore.errorOut(5, "-s must taked values of either 1 or 2");
			optParse(1);

		if args.branch_site not in [0,1,2]:
			gwctcore.errorOut(6, "-b must take values of 0, 1, or 2");
			optParse(1);

		if args.anc_opt not in [0,1]:
			gwctcore.errorOut(7, "-a must take values of 1 or 0");
			optParse(1);

		if args.verbosity not in [0,1]:
			gwctcore.errorOut(8, "-v must take values of either 1 or 0");
			optParse(1);

		if args.log_opt not in [0,1]:
			gwctcore.errorOut(9, "-l must take values of either 1 or 0");
			optParse(1);

		if args.logdir_suffix == None:
			args.logdir_suffix = "";

		return args.input, args.paml_path, args.tree_file, args.prune_opt, args.paml_seqtype, args.branch_site, args.anc_opt, args.verbosity, args.log_opt, args.logdir_suffix;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

#####

#def core.logCheck(lopt, lfilename, outline):
#	if lopt == 1:
#		core.printWrite(lfilename, outline);
#	else:
#		print outline;

############################################
#Main Block
############################################

ins, ppath, treefile, prune, seqtype, bsopt, aopt, v, l, outdir_suffix = optParse(0);

starttime = gwctcore.getLogTime();

if not os.path.isdir(ppath):
	gwctcore.errorOut(10, "-c must be a valid directory path");
	optParse(1);

if os.path.isfile(ins):
	fileflag = 1;
	indir = os.path.dirname(os.path.realpath(ins));
	indir, script_outdir = gwctcore.getOutdir(indir, "run_codeml", starttime, outdir_suffix);
	outdir = os.path.join(script_outdir, "codeml_out");
	if aopt == 1:
		ancdir = os.path.join(script_outdir, "anc_seqs_fa");
	filelist = [ins];

else:
	fileflag = 0;
	indir, script_outdir = gwctcore.getOutdir(ins, "run_codeml", starttime, outdir_suffix);
	outdir = os.path.join(script_outdir, "codeml_out");
	filelist = os.listdir(indir);
	if aopt == 1:
		ancdir = os.path.join(script_outdir, "anc_seqs_fa");

print gwctcore.getTime() + " | Creating main output directory:\t" + script_outdir;
os.system("mkdir " + script_outdir);

logfilename = os.path.join(script_outdir, "run_codeml.log");
logfile = open(logfilename, "w");
logfile.write("");
logfile.close();

gwctcore.logCheck(l, logfilename, "=======================================================================");
gwctcore.logCheck(l, logfilename, "\t\t\tRunning codeml");
gwctcore.logCheck(l, logfilename, "\t\t\t" + gwctcore.getDateTime());
if fileflag == 1:
	gwctcore.logCheck(l, logfilename, "INPUT    | Making tree from file:\t\t" + ins);
else:
	gwctcore.logCheck(l, logfilename, "INPUT    | Running codeml on all files in:\t" + indir);
gwctcore.logCheck(l, logfilename, "INFO     | PAML path set to:\t\t\t" + ppath);
if treefile != "":
	gwctcore.logCheck(l, logfilename, "INFO     | Using tree from file:\t\t" + treefile);
else:
	gwctcore.logCheck(l, logfilename, "INFO     | No tree file specified. codeml will infer a tree for each gene.");
if prune == 1:
	gwctcore.logCheck(l, logfilename, "INFO     | Pruning the tree for each gene.");
if bsopt == 0:
	gwctcore.logCheck(l, logfilename, "INFO     | Not doing branch-site test.");
if bsopt == 1:
	gwctcore.logCheck(l, logfilename, "INFO     | Doing NULL runs of branch-site test.");
if bsopt == 2:
	gwctcore.logCheck(l, logfilename, "INFO     | Doing ALT runs of branch-site test.");
if seqtype == 1:
	gwctcore.logCheck(l, logfilename, "INFO     | Seqtype set to codons. (1)");
if seqtype == 2:
	gwctcore.logCheck(l, logfilename, "INFO     | Seqtype set to amino acids. (2)");
if aopt == 0:
	gwctcore.logCheck(l, logfilename, "INFO     | Not performing ancestral reconstruction.");
elif aopt == 1:
	gwctcore.logCheck(l, logfilename, "INFO     | Saving ancestral reconstructions and probabilities.");
if v == 1:
	gwctcore.logCheck(l, logfilename, "INFO     | Printing all codeml output to the screen.");
else:
	gwctcore.logCheck(l, logfilename, "INFO     | Silent mode. Not printing codeml output to the screen.");
gwctcore.logCheck(l, logfilename, "OUTPUT   | An output directory has been created within the input directory called:\t" + script_outdir);
gwctcore.logCheck(l, logfilename, "-------------------------------------");
#sys.exit();
if not os.path.exists(outdir):
	gwctcore.logCheck(l, logfilename, gwctcore.getTime() + " | Creating codeml output directory:\t" + outdir);
	cmd = "mkdir " + outdir;
	os.system(cmd);

if aopt == 1:
	if not os.path.exists(ancdir):
		gwctcore.logCheck(l, logfilename, gwctcore.getTime() + " | Creating directory to pass ancestral sequences and trees:\t" + ancdir);
		cmd = "mkdir " + ancdir;
		os.system(cmd);

if prune == 1:
	gwctcore.logCheck(l, logfilename, gwctcore.getTime() + " | Retrieving tree info...");
	td, tree = gwctree.treeParse(open(treefile, "r").read().replace("\n",""),0);

	tips = [];
	for node in td:
		if td[node][2] == 'tip':
			tips.append(node);

gwctcore.logCheck(l, logfilename, gwctcore.getTime() + " | Starting codeml runs...\n");
if v == 0:
	codeml_logfile = os.path.join(script_outdir, "codeml.log");

ctlfilename = "codeml.ctl";

i = 0;
numbars = 0;
donepercent = [];

for each in filelist:
	if each.find(".fa") == -1:
		continue;

	if v == 0 and fileflag == 0:
		numbars, donepercent = gwctcore.loadingBar(i, len(filelist), donepercent, numbars);
	i = i + 1;

	if fileflag == 1:
		infilename = each;

	else:
		infilename = os.path.join(indir, each);
	gid = each[:each.index(".")];
	run_outdir = os.path.join(outdir, gid);
	if not os.path.exists(run_outdir):
		os.system("mkdir " + run_outdir);

	if prune == 1:
		seqs = gwctcore.fastaGetDict(infilename);
		to_prune = [];
		for tip in tips:
			if (">" + tip) not in seqs:
				to_prune.append(tip);

		nw_cmd = "nw_prune " + treefile + " ";
		for tip in to_prune:
			nw_cmd = nw_cmd + tip + " ";
		nw_cmd = nw_cmd + " > pruned.tre";
#		print nw_cmd;
		os.system(nw_cmd);

	ctlFile = open(ctlfilename, "w");

	inline = "seqfile = " + infilename + "\n";
	ctlFile.write(inline);
	if treefile != "":
		if prune == 1:
			treeline = "treefile = pruned.tre\n";
		else:
			treeline = "treefile = " + treefile + "\n";
		ctlFile.write(treeline);
	outline = "outfile = " + os.path.join(run_outdir, gid + ".out") + "\n\n";
	ctlFile.write(outline);

	ctlFile.write("noisy = 3\n");
	ctlFile.write("verbose = 0\n");
	if treefile != "":
		ctlFile.write("runmode = 0\n\n");
	else:
		ctlFile.write("runmode = 2\n\n");

	if seqtype == 1:
		ctlFile.write("seqtype = 1\n");
	elif seqtype == 2:
		ctlFile.write("seqtype = 2\n");

	ctlFile.write("CodonFreq = 2\n");
	ctlFile.write("clock = 0\n");
	ctlFile.write("aaDist = 0\n");
	ctlFile.write("aaRatefile = " + os.path.join(ppath, "dat", "wag.dat") + "\n");
	ctlFile.write("model = 2\n\n");

	if bsopt in [1,2]:
		ctlFile.write("NSsites = 2\n\n");
	else:
		ctlFile.write("NSsites = 0\n\n");

	ctlFile.write("icode = 0\n");
	ctlFile.write("fix_kappa = 0\n");
	ctlFile.write("kappa = 3\n");
	if bsopt == 1:
		ctlFile.write("fix_omega = 1\n");
	else:
		ctlFile.write("fix_omega = 0\n");
	ctlFile.write("omega = 1\n\n");

	ctlFile.write("fix_alpha = 1\n");
	ctlFile.write("alpha = 0\n");
	ctlFile.write("Malpha = 0\n");
	ctlFile.write("ncatG = 10\n\n");

	ctlFile.write("getSE = 0\n");
	if aopt == 1:
		ctlFile.write("RateAncestor = 1\n");
	ctlFile.write("Small_Diff = .5e-6\n");

	ctlFile.close();

	codeml_cmd = os.path.join(ppath, "bin", "codeml " + ctlfilename);
	if v == 0:
		if os.path.isfile(ins):
			codeml_cmd = codeml_cmd + " >> " + codeml_logfile;
		else:
			codeml_cmd = codeml_cmd + " >> " + codeml_logfile;
	if v == 1 or fileflag == 1:
		gwctcore.logCheck(l, logfilename, gwctcore.getTime() + " | codeml Call:\t" + codeml_cmd);
	else:
		lfile = open(logfilename, "a");
		lfile.write(gwctcore.getTime() + " | codeml Call for " + infilename + ":\t" + codeml_cmd + "\n");
		lfile.close();
	os.system(codeml_cmd);

	if aopt != 0:
		anc_seqfile = os.path.join(ancdir, gid + "_anc.fa");
		anc_probfile = os.path.join(ancdir, gid + "_ancprobs.fa");
		anc_treefile = os.path.join(ancdir, gid + "_anc.tre");
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
				if aopt == 1:
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
gwctcore.logCheck(l, logfilename, gwctcore.getTime() + " | Done!");
gwctcore.logCheck(l, logfilename, gwctcore.getLogTime());
gwctcore.logCheck(l, logfilename, "=======================================================================");
