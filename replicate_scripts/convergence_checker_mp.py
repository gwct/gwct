import sys, re, os, random, argparse
import multiprocessing as mp
import core

#MM spec IDs:
#"odoRosDiv1 lepWed1 turTru2 orcOrc1 triMan1"
#Initial background IDs:
#"hg19 rheMac3 calJac3 mm10 rn5 vicPac2 bosTau7 canFam3 loxAfr3 papHam1 monDom5"
#All specs:
#"Human":"hg19","Chimp":"panTro4","Gorilla":"gorGor3","Orangutan":"ponAbe2","Gibbon":"nomLeu3","Rhesus":"rheMac3","Crab-eating macaque":"macFas5","Baboon":"papHam1","Green monkey":"chlSab1","Marmoset":"calJac3","Squirrel monkey":"saiBol1","Bushbaby":"otoGar3","Chinese tree shrew":"tupChi1","Squirrel":"speTri2","Lesser Egyptian jerboa":"jacJac1","Prairie hamster":"micOch1","Chinese hamster":"criGri1","Golden hamster":"mesAur1","Mouse":"mm10","Rat":"rn5","Naked mole-rat":"hetGla2","Guinea pig":"cavPor3","Chinchilla":"chiLan1","Brush-tailed rat":"octDeg1","Rabbit":"oryCun2","Pika":"ochPri3","Pig":"susScr3","Alpaca":"vicPac2","Bactrian camel":"camFer1","Dolphin":"turTru2","Killer whale":"orcOrc1","Tibetan antelope":"panHod1","Cow":"bosTau7","Sheep":"oviAri3","Domestic goat":"capHir1","Horse":"equCab2","White rhinocerous":"cerSim1","Cat":"felCat5","Dog":"canFam3","Ferret":"musFur1","Panda":"ailMel1","Pacific walrus":"odoRosDiv1","Weddell seal":"lepWed1","Black flying-fox":"pteAle1","Megabat":"pteVam1","David's myotis bat":"myoDav1","Microbat":"myoLuc2","Big brown bat":"eptFus1","Hedgehog":"eriEur2","Shrew":"sorAra2","Star-nosed mole":"conCri1","Elephant":"loxAfr3","Cape elephant shrew":"eleEdw1","Manatee":"triMan1","Cape golden mole":"chrAsi1","Tenrec":"echTel2","Aardvark":"oryAfe1","Armadillo":"dasNov3","Opossum":"monDom5","Tasmanian devil":"sarHar1","Wallaby":"macEug2","Platypus":"ornAna1"

############################################
#Function Definitions
############################################
def IO_fileParse():
#This function handles the command line options.

	parser = argparse.ArgumentParser();

	parser.add_argument("-i", dest="input_directory", help="A directory containing FASTA formatted alignments you wish to check for convergence.");
	parser.add_argument("-t", dest="target_species", help="A space delimited list of species to look for convergent sites in. The entire list should be contained in quotes and the format of the species names must match that of the alignments.");
	parser.add_argument("-b", dest="background_species", help="A space delimited list of species to define as the background. The entire list should be contained in quotes and the format of the species names must match that of the alignments.");
	parser.add_argument("-r", dest="random_opt", help="If you wish to randomly select (background) species to check for convergent sites, enter the number of species you wish to select with this option.", type=int, default=0);
	parser.add_argument("-c", dest="num_reps", help="The number of times to perform the randomization test.", type=int, default=1);
	parser.add_argument("-d", dest="conv_divergence", help="Check for convergent divergent sites by setting this to 1.", type=int, default=0);
	parser.add_argument("-p", dest="num_threads", help="Multiple random replicates can be run in parallel. This sets the number of threads to be used. The quotient of -c and -p MUST be a whole number. Default: 1", type=int, default=1);
	parser.add_argument("-o", dest="output_file", help="Output file name for convergent genes/sites.");

	args = parser.parse_args();

	if args.input_directory == None or args.target_species == None or args.output_file == None:
		parser.print_help();
		sys.exit();

	if args.conv_divergence not in [0,1]:
		print " ------------------------------------------------";
		print "|**Error 1: -d must take values of either 0 or 1 |";
		print " ------------------------------------------------";
		parser.print_help();
		sys.exit();

	if args.random_opt != 0:
		if args.background_species != None:
			print " --------------------------------------------------------------------------------------------------------------------";
			print "|*Warning: With -r not set to 0, the background species will be chosen randomly. Your input species will be ignored. |"
			print " --------------------------------------------------------------------------------------------------------------------";
		bs = "";
		if args.num_reps == None:
			print " ------------------------------------------------------------";
			print "|**Error 2: With -r set to 1, a value for -c must be entered |";
			print " ------------------------------------------------------------";
			parser.print_help();
			sys.exit();
		elif (args.num_reps % args.num_threads) != 0:
			print " ---------------------------------------------------------------------------------------------";
			print "|**Error 3: When running random replicates, the quotient of -c and -p MUST be a whole number. |";
			print " ---------------------------------------------------------------------------------------------";
			parser.print_help();
			sys.exit();
		nt = args.num_threads;
	elif args.random_opt == 0:
		if args.background_species == 1:
			print " ----------------------------------------------------------------";
			print "|**Error 4: While -r is set to 0, -b must be defined by the user |";
			print " ----------------------------------------------------------------";
			parser.print_help();
			sys.exit();
		else:
			bs = args.background_species.split(" ");
		if args.num_threads > 1:
			nt = 1;
			print " ----------------------------------------------------------";
			print "|*Warning: Only 1 thread will be used while -r is set to 0. |"
			print " ----------------------------------------------------------";

	ts = args.target_species.split(" ");

	return args.input_directory, ts, bs, args.random_opt, args.num_reps, args.conv_divergence, nt, args.output_file;

#######

def remGapMiss(oldlist):
	newlist = [];
	for aa in oldlist:
		if aa in ["-", "X"]:
			continue;
		newlist.append(aa);
	return newlist;

#######
def convCheck(cur_c, c, ropt, targets, backgrounds, d, ins, outs):
#	cur_c = 0;
	init_c = cur_c+1;
	while cur_c < c:
		#if c > 1:
		if ropt != 0:
			outfilename = outs + "_" + str(cur_c+1) + ".txt";
		else:
			outfilename = outs + ".txt";

		if ropt != 0:
			#backgrounds = ["hg19", "rheMac3", "calJac3", "mm10", "rn5", "vicPac2", "bosTau7", "canFam3", "loxAfr3", "papHam1", "monDom5"];
			backgrounds = [];
			cur_r = len(backgrounds);
			while cur_r < ropt:
				chosenspec = random.choice(all_specs.values());

				if chosenspec not in targets and chosenspec not in backgrounds:
					backgrounds.append(chosenspec);
					cur_r = cur_r + 1;

		outfile = open(outfilename, "w");

		outfile.write("# ==============================================================================================\n");
		outfile.write("# \t\t\tConvergence testing\n");
		outfile.write("# \t\t\t" + core.getDateTime() + "\n");
		outfile.write("# Using alignments in:\t\t" + indir + "\n");
		outfile.write("# Target species:\t\t\t" + ", ".join(targets) + "\n");
		if ropt != 0:
			outfile.write("# Randomly choosing " + str(r) + " background species and performing " + str(c) + " replicate tests for convergence.\n");
			outfile.write("# This is replicate number " + str(cur_c+1) + "\n");
		outfile.write("# Background species:\t\t" + ", ".join(backgrounds) + "\n");
		outfile.write("# Writing output to:\t\t\t" + outfilename + "\n");
		if d == 0:
			outfile.write("# Checking for convergent sites.\n");
		elif d == 1:
			outfile.write("# Checking for divergent sites.\n");
		outfile.write("# ---------------------------------------------\n");
		#sys.exit();
		#cur_c = cur_c + 1;
		#continue;
		aligns = os.listdir(ins);

		numbars = 0;
		donepercent = [];
		count = len(aligns);
		i = 0;
		numsites = 0;
		totgenes = 0;
		outfile.write("# " + core.getTime() + " Starting Scan...\n");
		outfile.write("# Site#\tChromosome\tGeneID\tAlignLen\tPosition\tTargetAlleles\tBackgroundAlleles\n");
		for align in aligns:
			#numbars, donepercent = core.loadingBar(i, count, donepercent, numbars);
			i = i + 1;

			if align.find(".fa") == -1:
				continue;

			#if i > 25:
			#	break;

			infilename = ins + align;
			#print align;
			gid = "_".join(align.split("_")[:2]);
			chrome = align[align.find("chr"):align.find("chr")+4]

			inseqs = core.fastaGetDict(infilename);

			num_targets_present = 0;
			num_bg_present = 0;
			for title in inseqs:
				if any(t in title for t in targets):
					num_targets_present = num_targets_present + 1;
				if any(b in title for b in backgrounds):
					num_bg_present = num_bg_present + 1;

			if num_targets_present == len(targets) and num_bg_present == len(backgrounds):
				#print "The following gene has all target and background species and will be checked:\t\t" + gid;
				totgenes = totgenes + 1;

				seqlen = len(inseqs[inseqs.keys()[0]]);
				#print "Alignment length\t\t", seqlen;

				t_alleles = {};
				b_alleles = {};

				for x in xrange(len(inseqs[inseqs.keys()[0]])):
					for title in inseqs:
						for t in targets:
							if t in title:
								t_alleles[t] = inseqs[title][x];
						for b in backgrounds:
							if b in title:
								b_alleles[b] = inseqs[title][x];

					t_states = t_alleles.values();

					t_gap = t_states.count("-");
					t_missing = t_states.count("X");
					t_stop = t_states.count("*");

					b_states = b_alleles.values();

					b_gap = b_states.count("-");
					b_missing = b_states.count("X");
					b_stop = b_states.count("*");

					t_final = remGapMiss(t_states);
					b_final = remGapMiss(b_states);

					#print t_alleles;
					#print t_states;
					#print t_gap;
					#print t_missing;
					#print t_stop;
					#print t_final;

					#print b_alleles;
					#print b_states;
					#print b_gap;
					#print b_missing;
					#print b_stop;
					#print b_final;

					if t_final == [] or b_final == []:
						continue;

					if d == 0:
						if len(t_final) == len(targets) and len(b_final) == len(backgrounds) and t_final.count(t_final[0]) == len(t_final) and t_final[0] not in b_final:
							numsites = numsites + 1;
							#print core.getTime() + " Convergent site found!";
							#print "Filename:\t\t" + align;
							#print "Chromosome:\t\t" + chrome;
							#print "Gene ID:\t\t" + gid;
							#print "Alignment length\t", seqlen;
							#print "Target alleles:\t\t" + "".join(t_final);
							#print "Background alleles:\t" + "".join(b_final);
							#print "---------------";
							outline = str(numsites) + "\t" + chrome + "\t" + gid + "\t" + str(seqlen) + "\t" + str(x+1) + "\t" + "".join(t_final) + "\t" + "".join(b_final) + "\n";
							outfile.write(outline);

							#sys.exit();

					elif d == 1:
						if len(t_final) == len(targets) and len(b_final) == len(backgrounds) and t_final.count(t_final[0]) != len(t_final) and b_final.count(b_final[0]) == len(b_final):
							if not any(t in b_final for t in t_final):
								numsites = numsites + 1;
								#print "\nDivergent site found!";
								#print "Filename:\t\t" + align;
								#print "Chromosome:\t\t" + chrome;
								#print "Gene ID:\t\t" + gid;
								#print "Alignment length\t", seqlen;
								#print t_final;
								#print b_final;
								outline = str(numsites) + "\t" + chrome + "\t" + gid + "\t" + str(seqlen) + "\t" + str(x+1) + "\t" + "".join(t_final) + "\t" + "".join(b_final) + "\n";
								outfile.write(outline);

		#pstring = "100.0% complete.";
		#sys.stderr.write('\b' * len(pstring) + pstring);
		outfile.write("\n# " + core.getTime() + " Done!\n");
		outfile.write("# Total sites found: " + str(numsites) + "\n");
		outfile.write("# Total genes checked: " + str(totgenes) + "\n");
		outfile.write("# ==============================================================================================");
		cur_c = cur_c + 1;
	if ropt != 0:
		print core.getTime() + " Replicates", init_c, "to", c, "complete.";

############################################
#Main Block
############################################

indir, tgs, bgs, r, reps, cd, num_t, outfix = IO_fileParse();
all_specs = {"Human":"hg19","Chimp":"panTro4","Gorilla":"gorGor3","Orangutan":"ponAbe2","Gibbon":"nomLeu3","Rhesus":"rheMac3","Crab-eating macaque":"macFas5","Baboon":"papHam1","Green monkey":"chlSab1","Marmoset":"calJac3","Squirrel monkey":"saiBol1","Bushbaby":"otoGar3","Chinese tree shrew":"tupChi1","Squirrel":"speTri2","Lesser Egyptian jerboa":"jacJac1","Prairie hamster":"micOch1","Chinese hamster":"criGri1","Golden hamster":"mesAur1","Mouse":"mm10","Rat":"rn5","Naked mole-rat":"hetGla2","Guinea pig":"cavPor3","Chinchilla":"chiLan1","Brush-tailed rat":"octDeg1","Rabbit":"oryCun2","Pika":"ochPri3","Pig":"susScr3","Alpaca":"vicPac2","Bactrian camel":"camFer1","Dolphin":"turTru2","Killer whale":"orcOrc1","Tibetan antelope":"panHod1","Cow":"bosTau7","Sheep":"oviAri3","Domestic goat":"capHir1","Horse":"equCab2","White rhinocerous":"cerSim1","Cat":"felCat5","Dog":"canFam3","Ferret":"musFur1","Panda":"ailMel1","Pacific walrus":"odoRosDiv1","Weddell seal":"lepWed1","Black flying-fox":"pteAle1","Megabat":"pteVam1","David's myotis bat":"myoDav1","Microbat":"myoLuc2","Big brown bat":"eptFus1","Hedgehog":"eriEur2","Shrew":"sorAra2","Star-nosed mole":"conCri1","Elephant":"loxAfr3","Cape elephant shrew":"eleEdw1","Manatee":"triMan1","Cape golden mole":"chrAsi1","Tenrec":"echTel2","Aardvark":"oryAfe1","Opossum":"monDom5","Platypus":"ornAna1"};

if r > (len(all_specs) - len(tgs)):
	print " -------------------------------------------------------------------------------------------------------------------";
	print "|**Error 5: The number of species to be chosen randomly exceeds the total number of species available as background |";
	print " -------------------------------------------------------------------------------------------------------------------";
	sys.exit();

print "==============================================================================================";
print "\t\t\tConvergence testing";
print "\t\t\t" + core.getDateTime();
print "Using alignments in:\t\t" + indir
print "Target species:\t\t\t" + ", ".join(tgs);
if r != 0:
	print "Choosing " + str(r) + " random background species and performing " + str(reps) + " replicate tests.";
	print "The task will be split into " + str(num_t) + " processes and background species will be chosen as the processes are split.";
else:
	print "Background species:\t\t" + ", ".join(bgs);
	print "Writing output to:\t\t" + outfix + ".txt";
if cd == 0:
	print "Checking for convergent sites.";
elif cd == 1:
	print "Checking for divergent sites.";
print "---------------------------------------------";

if num_t > 1:
	bgs = [];
	processes = []
	output = mp.Queue()
	reps_per_t = reps / num_t;
	print core.getTime() + " Generating function calls...";
	print "Function calls and arguments are as follows:";
	print "convCheck([start replicate], [stop replicate], [# random species], [target species], [background species], [convergence/divergence option], [input directory], [output prefix])";
	print "----------";
	x = 0;
	while x < num_t:
		start_rep = reps_per_t * x;
		stop_rep = reps_per_t * (x+1);
		print core.getTime() + " Call " + str(x+1) + ":  convCheck(" + str(start_rep) + " , " + str(stop_rep) + " , " + str(r) + " , " + str(tgs) + " , " + str(bgs) + " , " + str(cd) + ", " + indir + " , " + outfix + ")";
		processes.append(mp.Process(target=convCheck,args=(start_rep,stop_rep,r,tgs,bgs,cd,indir,outfix)));
		x = x + 1;
	print "----------";
	##print processes;
	##convCheck(
	for p in processes:
		print core.getTime() + " start", p;
		p.start();
	for p in processes:
		print core.getTime() + " join", p;
		p.join();
	#print "Results:";
	#results = [output.get() for p in processes];
	#print results
	#sys.exit();

else:
	if r != 0:
		bgs = [];
		convCheck(0,1,r,tgs,bgs,cd,indir,outfix);
	else:
		convCheck(0,1,r,tgs,bgs,cd,indir,outfix);
print core.getTime() + " Done!";
print "==============================================================================================";
