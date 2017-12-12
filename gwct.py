#!/usr/bin/python
#############################################################################
#This is the IO wrapper of the GWCT. As input, you must have ancestral
#reconstructions from PAML as input. This would be the output directory of
#a run_codeml run (run_codeml is currently in the CORE repository).
#Gregg Thomas
#September 2015.
#############################################################################

import sys, os, platform, argparse, copy, math, lib.gwctcore as gc, lib.gwctree as gt, lib.convergence as convergence
import multiprocessing as mp

############################################
#Function Definitions
############################################
def optParse():
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Convergent substitution counting. Dependencies: ancestral sequences from PAML formatted with gwct_codeml.py");

	parser.add_argument("-i", dest="input", help="Input directory. This is the output directory of a gwct_codeml.py run and should include the folder 'anc-seqs-fa'. Will count convergent, divergent, and unique substitutions, or make pairwise comparisons if -w is set.");
	parser.add_argument("-t", dest="target_specs", help="A list contained in quotes. Spaces separate groups, commas separate species within groups. E.g: \"triMan1 lepWed1,odoRosDiv1 orcOrc1,turTru2\". If a target species is missing from a given alignment, it will be removed from the list of target species.", default="");
	parser.add_argument("-w", dest="pairwise_opt", help="Option to tell the program to simply do all pairwise comparisons of tip branches (1) or pairwise comparisons of all branches (2) and make a C/D graph. With -w set, unique substitutions will not be counted. Default: 0, do not do pairwise comparisons", type=int, default=0);
	parser.add_argument("-p", dest="prob_thresh", help="A probability threshold to only retrieve convergent sites with probabilities greater than or equal to. Set to 0 for no threshold. Default: 0", type=float, default=0);
	parser.add_argument("-n", dest="number_threads", help="The number of threads on which to run the job. NOTE: One thread will be reserved for the main process, so by entering '4' here, only 3 threads will be utilized on the data. Entering '2' is equivalent to running the data on 1 thread. Default: 1", type=int, default=1);
	parser.add_argument("-o", dest="output", help="Desired output directory. If none is entered, will be determined automatically.", default=False);

	args = parser.parse_args();

	if args.input == None or not os.path.isdir(args.input):
		sys.exit(gc.errorOut(1, "-i must be a valid directory path."));
	if args.pairwise_opt not in [0,1,2]:
		sys.exit(gc.errorOut(2, "-w must take values of 0, 1, or 2."));
	if args.target_specs == "" and args.pairwise_opt == 0:
		sys.exit(gc.errorOut(3, "With -w set to 0, -t must be defined."));	
	if args.target_specs != "" and args.pairwise_opt != 0:
		sys.exit(gc.errorOut(4, "Only one of -t and -w should be set."));
	if args.prob_thresh < 0.0 or args.prob_thresh > 1.0:
		sys.exit(gc.errorOut(5, "-p can only take values between 0.0 and 1.0."));
	if args.number_threads <= 0:
		sys.exit(gc.errorOut(8, "-n must be a positive, non-zero integer."));

	return args.input, args.target_specs, args.pairwise_opt, args.prob_thresh, args.number_threads, args.output;

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
	p = arglist[2];
	threads = arglist[3];
	results_dict = {};

	for filename in filelist_func:
		print filename;

		if "-ancprobs.fa" in filename:
			gid = filename[:filename.index("-ancprobs.fa")];
			treefilename = os.path.join(ancdir, gid + "-anc.tre");
		else:
			gid = filename[:filename.index("_ancprobs.fa")];
			treefilename = os.path.join(ancdir, gid + "_anc.tre");			
		infilename = os.path.join(ancdir, filename);
		
		tree = open(treefilename,"r").read().strip();
		tree_dict, new_tree, root = gt.treeParse(tree);
		#gene = "-".join(gid.split("_")[:2]);
		#chromosome = gid[gid.find("chr"):gid.find("chr")+4]
		
		if orig_targets != "":
			results_key = str(orig_targets);
			if results_key not in results_dict:
				results_dict[results_key] = [[],[],[]];
			targets = copy.deepcopy(orig_targets);
			#Resets the targets for each gene.

			results_dict = convergence.convCheck(infilename, results_dict, results_key, targets, prob_thresh, gid, tree_dict, pairwise);
			#Checking for convergent sites

		else:
			target_nodes = getTargs(tree_dict, p);
			for targets in target_nodes:
				if tree_dict[targets[0]][1] == targets[1] or tree_dict[targets[1]][1] == targets[0]:
					continue;
				# If one node is the ancestor of the other, skip this comparison.
				node_key = "";
				for n in targets:
					if "_" in n:
						node_key += n[n.index("_")+1:];
					else:
						node_key += tree_dict[n][3];
					if n == targets[0]:
						node_key += "-";

				if node_key not in results_dict:
					results_dict[node_key] = [[],[],[]];
	
				print targets;
				targets = [[t] for t in targets];
				print targets;

				results_dict = convergence.convCheck(infilename, results_dict, node_key, targets, prob_thresh, gid, tree_dict, pairwise);

	return results_dict;


############################################
#Main Block
############################################

starttime = gc.getLogTime();
indir, orig_targets, pairwise, prob_thresh, num_threads, output = optParse();

if platform.system() == "Windows" and num_threads != 1:
	print "\n** Warning! Multi-processing not currently supported on Windows. Switching to serial version (1 process).\n";
	num_threads = 1;

outdir = gc.defaultOut(output, indir, "-gwct");
print gc.getTime() + " | + Creating main output directory:\t" + outdir;
os.system("mkdir " + outdir);

ancdir = os.path.join(indir, "anc-seqs-fa");
if not os.path.isdir(ancdir):
	ancdir = os.path.join(indir, "anc_seqs_fa");
if not os.path.isdir(ancdir):
	sys.exit(gc.errorOut(9, "Cannot find anc-seqs-fa directory within input directory."));
filelist = os.listdir(ancdir);
filelist = [f for f in filelist if "ancprobs" in f];
num_files = len(filelist);
if prob_thresh == 0:
	convfilename = os.path.join(outdir, "conv-sites.txt");
	divfilename = os.path.join(outdir, "div-sites.txt");
elif prob_thresh != 0:
	convfilename = os.path.join(outdir, "conv-sites-" + str(prob_thresh) + ".txt");
	divfilename = os.path.join(outdir, "div-sites-" + str(prob_thresh) + ".txt");
uniqfilename = os.path.join(outdir, "unique-sites.txt");
main_header = "# GeneID\tAlignLen\tPosition\tAncAlleles\tTargetAlleles\n";
uniq_header = "# GeneID\tAlignLen\tPosition\tBackgroundAlleles\tTargetAlleles\n";
if pairwise == 0:
	gc.filePrep(convfilename, main_header);
	gc.filePrep(divfilename, main_header);
	gc.filePrep(uniqfilename, uniq_header);

logfilename = os.path.join(outdir, "gwct.log");
gc.filePrep(logfilename, "");

sp = 65;
gc.printWrite(logfilename, "# =======================================================================");
gc.printWrite(logfilename, "#\t\tGWCT: Genome-Wide Convergence Tester");
gc.printWrite(logfilename, "#\t\t\t" + gc.getDateTime());
gc.printWrite(logfilename, gc.spacedOut("# INPUT  | Using ancestral reconstructions in:", sp) + indir);
if orig_targets != "":
	orig_targets = [t.split(",") for t in orig_targets.split(" ")];
	gc.printWrite(logfilename, gc.spacedOut("# INFO   | Target species:", sp) + str(orig_targets));
elif pairwise == 1:
	pw_dir = outdir + "pairwise-tips/"
	pw_dir = os.path.join(outdir, "pairwise-tips");
	gc.printWrite(logfilename, "# INFO   | Performing pairwise comparisons of TIP branches");
elif pairwise == 2:
	pw_dir = os.path.join(outdir, "pairwise-all");
	gc.printWrite(logfilename, "# INFO   | Performing pairwise comparisons of ALL branches");
if prob_thresh > 0:
	gc.printWrite(logfilename, gc.spacedOut("# INFO   | Ancestral probability threshold:", sp) + str(prob_thresh));
gc.printWrite(logfilename, gc.spacedOut("# INFO   | Number of threads:", sp) + str(num_threads));
gc.printWrite(logfilename, gc.spacedOut("# OUTPUT | Output directory created within input directory:", sp) + outdir);
gc.printWrite(logfilename, gc.spacedOut("# INFO   | Checking for convergent, divergent, and unique sites", sp));
gc.printWrite(logfilename, "# ---------------------------------------------\n");

if pairwise != 0:
	print gc.getTime() +  " | + Creating pairwise output directory.\n";
	os.system("mkdir " + pw_dir);
	target_counts = {};

################################

###SPLIT UP THE GENES, MAKE THE ARGS, CALL THE FUNCTION
if num_threads == 1:
	main_args = [filelist, orig_targets, pairwise, num_threads]
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
		i += 1;
		if i > files_per_split:
			i = 1;
			split_num = split_num + 1;

	pool = mp.Pool(processes=(num_threads-1));

	args = []
	for files in filesplit:
		args.append([files, orig_targets, pairwise, num_threads]);

	results = pool.map(splitThreads,args)

	print "Results:";

	main_results_dict = {};
	for proc in results:
		for key in proc:
			if key not in main_results_dict:
				main_results_dict[key] = [[],[],[]];
			main_results_dict[key][0] = main_results_dict[key][0] + proc[key][0];
			main_results_dict[key][1] = main_results_dict[key][1] + proc[key][1];
			if pairwise == 0:
				main_results_dict[key][2] = main_results_dict[key][2] + proc[key][2];

##################

if orig_targets != "":
	main_key = str(orig_targets);
	open(convfilename,"a").writelines(main_results_dict[main_key][0]);
	conv_sites = len(main_results_dict[main_key][0]);
	open(divfilename,"a").writelines(main_results_dict[main_key][1]);
	div_sites = len(main_results_dict[main_key][1]);
	open(uniqfilename,"a").writelines(main_results_dict[main_key][2]);
	uniq_sites = len(main_results_dict[main_key][2]);

else:
	conv_sites, div_sites, uniq_sites = 0,0,0;
	conv_sites_vect, div_sites_vect = [],[];

	lfile = open(logfilename, "a");
	pairwise_header = "Node#1\tNode#2\t#Convergent\t#Divergent\t#Unique\n";
	lfile.write(pairwise_header);
	for pair in main_results_dict:
		pairconvfile = os.path.join(pw_dir, pair + "-conv.txt");
		pairdivfile = os.path.join(pw_dir, pair + "-div.txt");
		
		gc.filePrep(pairconvfile, main_header);
		gc.filePrep(pairdivfile, main_header);
		
		open(pairconvfile, "a").writelines(main_results_dict[pair][0]);
		open(pairdivfile, "a").writelines(main_results_dict[pair][1]);
		
		cur_conv_num = len(main_results_dict[pair][0]);
		cur_div_num = len(main_results_dict[pair][1]);
		
		outline = pair + "\t" + str(cur_conv_num) + "\t" + str(cur_div_num);
		if pairwise in [0,1]:
			pairuniqfile = os.path.join(pw_dir, pair + "-uniq.txt");
			gc.filePrep(pairuniqfile, uniq_header);
			open(pairuniqfile, "a").writelines(main_results_dict[pair][2]);
			cur_uniq_num = len(main_results_dict[pair][2]);
			outline += "\t" + str(cur_uniq_num);
			uniq_sites += cur_uniq_num;
		lfile.write(outline + "\n");
		conv_sites_vect.append(cur_conv_num);
		div_sites_vect.append(cur_div_num);
		conv_sites += cur_conv_num;
		div_sites += cur_div_num;
	lfile.close();

	open(os.path.join(outdir, "crtmp.txt"),"w").write(" ".join(map(str, conv_sites_vect)));
	open(os.path.join(outdir, "drtmp.txt"),"w").write(" ".join(map(str, div_sites_vect)));
	gc.printWrite(logfilename, "\n\n# " + gc.getTime() + " | Calling R to plot C/D ratios\n");
	try:
		print "-Begin R screen output-";
		r_cmd = "rscript " + os.path.join(os.path.dirname(os.path.realpath(__file__)), "lib", "cd_ratio_reg.r") + " " + outdir + " crtmp.txt drtmp.txt";
		os.system(r_cmd);
		print "--End R screen output--";
	except:
		print gc.printWrite(logfilename, "** ERROR RUNNING RSCRIPT... SKIPPING C/D PLOT GENERATION.");
	gc.printWrite(logfilename, "\n# " + gc.getTime() + " | Removing temporary files");
	os.system("rm " + os.path.join(outdir, "crtmp.txt"));
	os.system("rm " + os.path.join(outdir, "drtmp.txt"));

gc.printWrite(logfilename, "\n# ---------------------------------------------");
gc.printWrite(logfilename, "\n# " + gc.getTime() + " Done!");
if type(conv_sites) == list:
	conv_sites = sum(conv_sites);
	div_sites = sum(div_sites);
gc.printWrite(logfilename, "# Total convergent sites found: " + str(conv_sites));
gc.printWrite(logfilename, "# Total divergent sites found: " + str(div_sites));
if pairwise in [0,1]:
	gc.printWrite(logfilename, "# Total unique sites found: " + str(uniq_sites));
gc.printWrite(logfilename, "# Total genes checked: " + str(num_files));
gc.printWrite(logfilename, "# ==============================================================================================");
