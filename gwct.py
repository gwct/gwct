#!/usr/bin/python
#############################################################################
#This is the IO wrapper of the GWCT. As input, you must have ancestral
#reconstructions from PAML as input. This would be the output directory of
#a run_codeml run (run_codeml is currently in the CORE repository).
#Gregg Thomas
#September 2015.
#############################################################################

import sys, os, argparse, copy, math, lib.gwctcore as gwctcore, lib.gwctree as gwctree, lib.convergence as convergence
import multiprocessing as mp

############################################
#Function Definitions
############################################
def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Convergent substitution counting. Dependencies: ancestral sequences from PAML");

	parser.add_argument("-i", dest="input", help="Input directory containing a single run_codeml output directory or codeml_combined directory.");
	parser.add_argument("-t", dest="target_specs", help="A list contained in quotes. Spaces separate groups, commas separate species within groups. E.g: \"triMan1 lepWed1,odoRosDiv1 orcOrc1,turTru2\"", default="");
	parser.add_argument("-u", dest="uniq_subs", help="A boolean to output the sites that are unique substitutions in the tip branches of interest (1) or not (0). Default: 0", type=int, default = 0);
	parser.add_argument("-w", dest="pairwise_opt", help="Option to tell the program to simply do all pairwise comparisons of tip branches (1) or pairwise comparisons of all branches (2) and make a C/D graph. Default: 0, do not do pairwise comparisons", type=int, default=0);
	parser.add_argument("-p", dest="prob_thresh", help="A probability threshold to only retrieve convergent sites with probabilities greater than or equal to. Set to 0 for no threshold. Default: 0", type=float, default=0);
	parser.add_argument("-n", dest="number_threads", help="The number of threads on which to run the job. NOTE: One thread will be reserved for the main process, so by entering '4' here, only 3 threads will be utilized on the data. Entering '2' is equivalent to running the data on 1 thread. Default: 1", type=int, default=1);
	parser.add_argument("-o", dest="output_suffix", help="The SUFFIX of the directory name to be created by the script. Something descriptive for your run.", default="");

	args = parser.parse_args();

	if errorflag == 0:
		if args.input == None:
			gwctcore.errorOut(1, "-i must be defined");
			optParse(1);

		if args.pairwise_opt not in [0,1,2]:
			gwctcore.errorOut(2, "-w must take values of 0, 1, or 2");
			optParse(1);

		if args.target_specs == "" and args.pairwise_opt == 0:
			gwctcore.errorOut(3, "With -w set to 0, -t must be defined");
			optParse(1);

		if args.target_specs != "" and args.pairwise_opt != 0:
			gwctcore.errorOut(4, "Only one of -t and -w should be set");
			optParse(1);

		if args.prob_thresh < 0.0 or args.prob_thresh > 1.0:
			gwctcore.errorOut(5, "-p can only take values between 0.0 and 1.0");
			optParse(1);

		if args.uniq_subs not in [0,1]:
			gwctcore.errorOut(6, "-u can only take values of 0 or 1");
			optParse(1);

		if args.uniq_subs == 1 and args.target_specs == "":
			gwctcore.errorOut(7, "-u can only be set when a set of target species is defined with -t");
			optParse(1);

		if args.uniq_subs == 1 and args.target_specs == None:
			gwctcore.errorOut(8, "With -u set to 1, -t must also be defined and -w must be 0");
			optParse(1);

		if args.number_threads <= 0:
			gwctcore.errorOut(9, "-n must be a positive, non-zero integer");
			optParse(1);

		if args.output_suffix == None:
			args.output_suffix = "";

		return args.input, args.target_specs, args.uniq_subs, args.pairwise_opt, args.prob_thresh, args.number_threads, args.output_suffix;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

##########
def getTargs(td, t_opt):
	target_list = []
	for t1 in td:
		if (t_opt == 1 and td[t1][2] != 'tip') or (t_opt == 2 and td[t1][2] == 'root'):
			continue;
		for t2 in td:
			if (t_opt == 1 and td[t2][2] != 'tip') or (t_opt == 2 and td[t2][2] == 'root') or t1 == t2:
				continue;
			if [t1,t2] not in target_list and [t2,t1] not in target_list:
				target_list.append([t1,t2]);

	return target_list;

############################################
#Threaded Block
############################################
#def splitThreads(filelist_func, orig_targets_func, u, p, threads):
def splitThreads(arglist):
	filelist_func = arglist[0];
	orig_targets = arglist[1];
	u = arglist[2];
	p = arglist[3];
	threads = arglist[4];
	results_dict = {};

	for filename in filelist_func:
		print filename;

		gid = filename[:filename.index("_ancprobs.fa")];
		gene = "_".join(gid.split("_")[:2]);
		chromosome = gid[gid.find("chr"):gid.find("chr")+4]
		infilename = os.path.join(indir, filename);
		treefilename = os.path.join(indir, gid + "_anc.tre");
		tree = open(treefilename,"r").read().replace("\n","");
		tree_dict, new_tree = gwctree.treeParse(tree);

		if orig_targets != "":
			results_key = str(orig_targets_func);
			if results_key not in results_dict:
				results_dict[results_key] = [[],[],[]];
			targets = copy.deepcopy(orig_targets_func);
			#Resets the targets for each gene.
			results_dict = convergence.convCheck(infilename, results_dict, results_key, targets, prob_thresh, chromosome, gene, tree_dict, u);
			#Checking for convergent sites

		else:
			target_nodes = getTargs(tree_dict, p);
			for targets in target_nodes:
				if tree_dict[targets[0]][1] == targets[1] or tree_dict[targets[1]][1] == targets[0]:
					continue;
				#If one node is the ancestor of the other, skip this comparison.

				node_key = "";
				for n in targets:
					if "_" in n:
						node_key = node_key + n[n.index("_")+1:];
					else:
						node_key = node_key + n;
					if n == targets[0]:
						node_key = node_key + "-";

				if node_key not in results_dict:
					results_dict[node_key] = [[],[],[]];

				results_dict = convergence.convCheck(infilename, results_dict, node_key, targets, prob_thresh, chromosome, gene, tree_dict, u);

	return results_dict;


############################################
#Main Block
############################################

starttime = gwctcore.getLogTime();
indir, orig_targets, unique, pairwise, prob_thresh, num_threads, suffix = optParse(0);

if suffix != "":
	suffix = "-" + suffix;
if sys.platform.find("win") != -1 and sys.platform != "darwin":
	rm_cmd = "del ";
else:
	rm_cmd = "rm ";

indir, script_outdir = gwctcore.getOutdir(indir, "gwct", starttime, suffix);
print indir;
inslist = os.listdir(indir);
for each in inslist:
	if each.find("run_codeml") != -1 or each.find("codeml_combined") != -1:
		#indir = indir + each + "/anc_seqs_fa/";
		indir = os.path.join(indir, each, "anc_seqs_fa");
		orig_filelist = os.listdir(indir);
		filelist = [f for f in orig_filelist if f.find("ancprobs") != -1];
		num_files = len(filelist);
		break;

print gwctcore.getTime() + " | Creating main output directory:\t" + script_outdir;
os.system("mkdir " + script_outdir);

logfilename = os.path.join(script_outdir, "gwct.log");
logfile = open(logfilename, "w");
logfile.write("");
logfile.close();
main_header = "# Chromosome\tGeneID\tAlignLen\tPosition\tAncAlleles\tTargetAlleles\n";
uniq_header = "# Chromosome\tGeneID\tAlignLen\tPosition\tBackgroundAlleles\tTargetAlleles\n";
l = 1;
sp = 65;

gwctcore.logCheck(l, logfilename, "# =======================================================================");
gwctcore.logCheck(l, logfilename, "#\t\tGWCT: Genome-Wide Convergence Tester");
gwctcore.logCheck(l, logfilename, "#\t\t\t" + gwctcore.getDateTime());
gwctcore.logCheck(l, logfilename, gwctcore.spacedOut("# INPUT  | Using ancestral reconstructions in:", sp) + indir);
if orig_targets != "":
	orig_targets = orig_targets.split(" ");
	for t in xrange(len(orig_targets)):
		if orig_targets[t].find(",") != -1:
			orig_targets[t] = orig_targets[t].split(",");
	if prob_thresh == 0:
		convfilename = os.path.join(script_outdir, "conv_sites.txt");
		divfilename = os.path.join(script_outdir, "div_sites.txt");
	elif prob_thresh != 0:
		convfilename = os.path.join(script_outdir, "conv_sites_" + str(prob_thresh) + ".txt");
		divfilename = os.path.join(script_outdir, "div_sites_" + str(prob_thresh) + ".txt");
	gwctcore.filePrep(convfilename, main_header);
	gwctcore.filePrep(divfilename, main_header);
	gwctcore.logCheck(l, logfilename, gwctcore.spacedOut("# INFO   | Target species:", sp) + str(orig_targets));
	if unique == 1:
		uniqfilename = os.path.join(script_outdir, "unique_sites.txt");
		gwctcore.filePrep(uniqfilename, uniq_header);
		gwctcore.logCheck(l, logfilename, "# INFO   | Counting unique substitutions on target branches");
	else:
		uniqfilename = "";
#	gwctcore.logCheck(l, logfilename, "# INFO   | Not performing pairwise comparisons");
elif pairwise == 1:
	pw_dir = script_outdir + "pairwise_tips/"
	pw_dir = os.path.join(script_outdir, "pairwise_tips");
	gwctcore.logCheck(l, logfilename, "# INFO   | Performing pairwise comparisons of TIP branches");
elif pairwise == 2:
	pw_dir = os.path.join(script_outdir, "pairwise_all");
	gwctcore.logCheck(l, logfilename, "# INFO   | Performing pairwise comparisons of ALL branches");
if prob_thresh > 0:
	gwctcore.logCheck(l, logfilename, gwctcore.spacedOut("# INFO   | Ancestral probability threshold:", sp) + str(prob_thresh));
gwctcore.logCheck(l, logfilename, gwctcore.spacedOut("# INFO   | Number of threads:", sp) + str(num_threads));
gwctcore.logCheck(l, logfilename, gwctcore.spacedOut("# OUTPUT | Output directory created within input directory:", sp) + script_outdir);
gwctcore.logCheck(l, logfilename, gwctcore.spacedOut("# INFO   | Checking for convergent and divergent sites", sp));
gwctcore.logCheck(l, logfilename, "# ---------------------------------------------\n");
#sys.exit()
if pairwise != 0:
	uniqfilename = "";
	print gwctcore.getTime() +  " | +Creating pairwise output directory.\n";
	os.system("mkdir " + pw_dir);
	target_counts = {};

################################

###SPLIT UP THE GENES, MAKE THE ARGS, CALL THE FUNCTION
if num_threads == 1:
	main_args = [filelist, orig_targets, unique, pairwise, num_threads]
	main_results_dict = splitThreads(main_args);
else:
	filesplit = [];
	files_per_split = math.ceil(float(num_files) / float(num_threads-1));
	split_num = 0;
	i = 1;
	for each in filelist:
		if i == 1:
			filesplit.append([]);
		filesplit[split_num].append(each);
		i = i + 1;
		if i > files_per_split:
			i = 1;
			split_num = split_num + 1;


	pool = mp.Pool(processes=(num_threads-1));

	args = []
	for files in filesplit:
		args.append([files, orig_targets, unique, pairwise, num_threads]);
	results = pool.map(splitThreads,args)

	print "Results:";

	main_results_dict = {};
	for proc in results:
		for key in proc:
			if key not in main_results_dict:
				main_results_dict[key] = [[],[],[]];
			main_results_dict[key][0] = main_results_dict[key][0] + proc[key][0];
			main_results_dict[key][1] = main_results_dict[key][1] + proc[key][1];
			if unique == 1:
				main_results_dict[key][2] = main_results_dict[key][2] + proc[key][2];

	#print main_results_dict;

if orig_targets != "":
	main_key = str(orig_targets);
	open(convfilename,"a").writelines(main_results_dict[main_key][0]);
	conv_sites = len(main_results_dict[main_key][0]);
	open(divfilename,"a").writelines(main_results_dict[main_key][1]);
	div_sites = len(main_results_dict[main_key][1]);
	if unique == 1:
		open(uniqfilename,"a").writelines(main_results_dict[main_key][2]);
		uniq_sites = len(main_results_dict[main_key][2]);


else:
	conv_sites = 0;
	div_sites = 0;
	conv_sites_vect = [];
	div_sites_vect = [];

	lfile = open(logfilename, "a");
	pairwise_header = "Node #1\tNode#2\t# Convergent\t#Divergent\n";
	lfile.write(pairwise_header);
	for pair in main_results_dict:
		pairconvfile = os.path.join(pw_dir, pair + "_conv.txt");
		pairdivfile = os.path.join(pw_dir, pair + "_div.txt");
		gwctcore.filePrep(pairconvfile, main_header);
		gwctcore.filePrep(pairdivfile, main_header);
		open(pairconvfile, "a").writelines(main_results_dict[pair][0]);
		open(pairdivfile, "a").writelines(main_results_dict[pair][1]);
		cur_conv_num = len(main_results_dict[pair][0]);
		cur_div_num = len(main_results_dict[pair][1]);
		outline = pair + "\t" + str(cur_conv_num) + "\t" + str(cur_div_num) + "\n";
		lfile.write(outline);
		conv_sites_vect.append(cur_conv_num);
		div_sites_vect.append(cur_div_num);
		conv_sites = conv_sites + cur_conv_num;
		div_sites = div_sites + cur_div_num;
	lfile.close();

	open(os.path.join(script_outdir, "crtmp.txt"),"w").write(" ".join(map(str, conv_sites_vect)));
	open(os.path.join(script_outdir, "drtmp.txt"),"w").write(" ".join(map(str, div_sites_vect)));
	gwctcore.logCheck(l, logfilename, "\n\n# " + gwctcore.getTime() + " | Calling R to plot C/D ratios\n");
	print "-Begin R screen output-";
	r_cmd = "rscript " + os.path.join(os.path.dirname(os.path.realpath(__file__)), "lib", "cd_ratio_reg.r") + " " + script_outdir + " crtmp.txt drtmp.txt";
	os.system(r_cmd);
	print "--End R screen output--";
	gwctcore.logCheck(l, logfilename, "\n# " + gwctcore.getTime() + " | Removing temporary files");
	os.system(rm_cmd + os.path.join(script_outdir, "crtmp.txt"));
	os.system(rm_cmd + os.path.join(script_outdir, "drtmp.txt"));


gwctcore.logCheck(l, logfilename, "\n# ---------------------------------------------");
gwctcore.logCheck(l, logfilename, "\n# " + gwctcore.getTime() + " Done!");
if type(conv_sites) == list:
	conv_sites = sum(conv_sites);
	div_sites = sum(div_sites);
gwctcore.logCheck(l, logfilename, "# Total convergent sites found: " + str(conv_sites));
gwctcore.logCheck(l, logfilename, "# Total divergent sites found: " + str(div_sites));
if unique == 1:
	gwctcore.logCheck(l, logfilename, "# Total unique sites found: " + str(uniq_sites));
gwctcore.logCheck(l, logfilename, "# Total genes checked: " + str(num_files));
gwctcore.logCheck(l, logfilename, "# ==============================================================================================");
