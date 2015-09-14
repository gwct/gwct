#############################################################################
#This is the IO wrapper of the GWCT. As input, you must have ancestral 
#reconstructions from PAML as input. This would be the output directory of
#a run_codeml run (run_codeml is currently in the CORE repository).
#Gregg Thomas
#September 2015.
#############################################################################

import sys, os, argparse, copy
sys.path.append(sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/lib-gwct/"))
import gwctlib, gwctree
import convergence

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
	#parser.add_argument("-o", dest="output", help="A name for an output file. The file will be placed in the run_codeml directory within the input directory.");

	args = parser.parse_args();

	if errorflag == 0:
		if args.input == None:
			gwctlib.errorOut(1, "-i must be defined");
			optParse(1);

		if args.pairwise_opt not in [0,1,2]:
			gwctlib.errorOut(2, "-w must take values of 0, 1, or 2");
			optParse(1);

		if args.target_specs == "" and args.pairwise_opt == 0:
			gwctlib.errorOut(3, "With -w set to 0, -t must be defined");
			optParse(1);

		if args.target_specs != "" and args.pairwise_opt != 0:
			gwctlib.errorOut(4, "Only one of -t and -w should be set");
			optParse(1);

		if args.prob_thresh < 0.0 or args.prob_thresh > 1.0:
			gwctlib.errorOut(5, "-p can only take values between 0.0 and 1.0");
			optParse(1);

		if args.uniq_subs not in [0,1]:
			gwctlib.errorOut(6, "-u can only take values of 0 or 1");
			optParse(1);

		if args.uniq_subs == 1 and args.target_specs == None:
			gwctlib.errorOut(7, "With -u set to 1, -t must also be defined and -w must be 0");
			optParse(1);
		
		return args.input, args.target_specs, args.uniq_subs, args.pairwise_opt, args.prob_thresh;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

##########
def getTargs(td, t_opt):
	target_list = []
	for t1 in tree_dict:
		if (t_opt == 1 and tree_dict[t1][9] != 'tip') or (t_opt == 2 and tree_dict[t1][9] == 'root'):
			continue;
		for t2 in tree_dict:
			if (t_opt == 1 and tree_dict[t2][9] != 'tip') or (t_opt == 2 and tree_dict[t2][9] == 'root') or t1 == t2:
				continue;
			if t1+"-"+t2 not in target_list and t2+"-"+t1 not in target_list:
				target_list.append([t1,t2]);

	return target_list;

############################################
#Main Block
############################################

starttime = gwctlib.getLogTime();
indir, orig_targets, unique, pairwise, prob_thresh = optParse(0);

indir, script_outdir = gwctlib.getOutdir(indir, "gt_conv", starttime);
inslist = os.listdir(indir);
for each in inslist:
	if each.find("run_codeml") != -1 or each.find("codeml_combined") != -1:
		indir = indir + each + "/anc_seqs_fa/";
		filelist = os.listdir(indir);
		break;

print gwctlib.getTime() + " | Creating main output directory:\t" + script_outdir;
os.system("mkdir " + script_outdir);

logfilename = script_outdir + "gt_conv.log";
logfile = open(logfilename, "w");
logfile.write("");
logfile.close();
main_header = "# Site#\tChromosome\tGeneID\tAlignLen\tPosition\tAncAlleles\tTargetAlleles\n";
l = 1;
sp = 65;

gwctlib.logCheck(l, logfilename, "=======================================================================");
gwctlib.logCheck(l, logfilename, "\t\t\tCounting convergent sites");
gwctlib.logCheck(l, logfilename, "\t\t\t" + gwctlib.getDateTime());
gwctlib.logCheck(l, logfilename, gwctlib.spacedOut("# INPUT  | Using ancestral reconstructions in:", sp) + indir);
if orig_targets != "":
	orig_targets = orig_targets.split(" ");
	for t in xrange(len(orig_targets)):
		if orig_targets[t].find(",") != -1:
			orig_targets[t] = orig_targets[t].split(",");
	if prob_thresh == 0:
		convfilename = script_outdir + "conv_sites.txt";
		divfilename = script_outdir + "div_sites.txt";
	elif prob_thresh != 0:
		convfilename = script_outdir + "conv_sites_" + str(prob_thresh) + ".txt";
		divfilename = script_outdir + "div_sites_" + str(prob_thresh) + ".txt";
	gwctlib.filePrep(convfilename, main_header);
	gwctlib.filePrep(divfilename, main_header);
	gwctlib.logCheck(l, logfilename, gwctlib.spacedOut("# INFO   | Target species:", sp) + str(orig_targets));
	if unique == 1:
		if prob_thresh == 0:
			convfilename = script_outdir + "unique_sites.txt";
		elif prob_thresh != 0:
			convfilename = script_outdir + "unique_sites_" + str(prob_thresh) + ".txt";
		gwctlib.filePrep(uniqfilename, "");
		gwctlib.logCheck(l, logfilename, "# INFO   | Counting unique substitutions on target branches");
#	gwctlib.logCheck(l, logfilename, "# INFO   | Not performing pairwise comparisons");
elif pairwise == 1:
	pw_dir = script_outdir + "pairwise_tips/"
	gwctlib.logCheck(l, logfilename, "# INFO   | Performing pairwise comparisons of TIP branches");
elif pairwise == 2:
	pw_dir = script_outdir + "pairwise_all/";
	gwctlib.logCheck(l, logfilename, "# INFO   | Performing pairwise comparisons of ALL branches");
if prob_thresh > 0:
	gwctlib.logCheck(l, logfilename, gwctlib.spacedOut("# INFO   | Ancestral probability threshold:", sp) + str(prob_thresh));
gwctlib.logCheck(l, logfilename, gwctlib.spacedOut("# OUTPUT | Output directory created within input directory:", sp) + script_outdir);
gwctlib.logCheck(l, logfilename, gwctlib.spacedOut("# INFO   | Checking for convergent and divergent sites.", sp));
gwctlib.logCheck(l, logfilename, "# ---------------------------------------------\n");

if pairwise != 0:
	print gwctlib.getTime() +  " | +Creating pairwise output directory.\n";
	os.system("mkdir " + pw_dir);
	target_counts = {};



################################

conv_sites = 0;
div_sites = 0;

numbars = 0;
donepercent = [];
numfiles = len(filelist);
i = 0;
j = 0;
for filename in filelist:
	numbars, percentdone = gwctlib.loadingBar(j, numfiles, donepercent, numbars);
	j = j + 1;
	if filename.find("ancprobs.fa") == -1:
		continue;
	i = i + 1;
	#For every gene, however if its not an ancestral probabilities file, skip it. j is total number of files, i is total number of genes.
	#print filename;

	gid = filename[:filename.index("_ancprobs.fa")];
	gene = "_".join(gid.split("_")[:2]);
	chromosome = gid[gid.find("chr"):gid.find("chr")+4]
	infilename = indir + filename;
	treefilename = indir + gid + "_anc.tre";
	tree = open(treefilename,"r").read().replace("\n","");
	tree_dict, new_tree = gwctree.treeParse(tree,2);
	#print new_tree;

	if orig_targets != "":
		targets = copy.deepcopy(orig_targets);
		#Resets the targets for each gene.
		conv_sites, div_sites = convergence.convCheck(infilename, convfilename, divfilename, targets, conv_sites, div_sites, prob_thresh, chromosome, gene, tree_dict);
		#Checking for convergent sites
	
	else:
		target_nodes = getTargs(tree_dict, pairwise);
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

			if node_key not in target_counts:
				target_counts[node_key] = [0,0];

			conv_sites = 0;
			div_sites = 0;

			pairconvfile = pw_dir + node_key + "_conv.txt";
			pairdivfile = pw_dir + node_key + "_div.txt";
			gwctlib.filePrep(pairconvfile, main_header);
			gwctlib.filePrep(pairdivfile, main_header);

			conv_sites, div_sites = convergence.convCheck(infilename, pairconvfile, pairdivfile, targets, conv_sites, div_sites, prob_thresh, chromosome, gene, tree_dict);
			target_counts[node_key][0] = target_counts[node_key][0] + conv_sites;
			target_counts[node_key][1] = target_counts[node_key][1] + div_sites;

	if i > 0:
		break;	

pstring = "100.0% complete.";
sys.stderr.write('\b' * len(pstring) + pstring);

if pairwise != 0:
	conv_sites = [];
	div_sites = [];
	lfile = open(logfilename, "a");
	pairwise_header = "Node #1\tNode#2\t# Convergent\t#Divergent\n";
	lfile.write(pairwise_header);
	for pair in target_counts:
		nodes = pair.split("-");
		outline = nodes[0] + "\t" + nodes[1] + "\t" + str(target_counts[pair][0]) + "\t" + str(target_counts[pair][1]) + "\n";
		conv_sites.append(target_counts[pair][0]);
		div_sites.append(target_counts[pair][1])
		lfile.write(outline);
	lfile.close();

gwctlib.logCheck(l, logfilename, "\n# ---------------------------------------------");
gwctlib.logCheck(l, logfilename, "\n# " + gwctlib.getTime() + " Done!");
if type(conv_sites) == list:
	conv_sites = sum(conv_sites);
	div_sites = sum(div_sites);
gwctlib.logCheck(l, logfilename, "# Total convergent sites found: " + str(conv_sites));
gwctlib.logCheck(l, logfilename, "# Total divergent sites found: " + str(div_sites));
gwctlib.logCheck(l, logfilename, "# Total genes checked: " + str(i));
gwctlib.logCheck(l, logfilename, "# ==============================================================================================");

for t in target_counts:
	print t, ":", target_counts[t];
















