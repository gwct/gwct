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
def optParse(errorflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser();

	parser.add_argument("-i", dest="input_directory", help="A directory containing FASTA formatted alignments you wish to check for convergence.");
	parser.add_argument("-s", dest="num_spec", help="The number of species to be chosen from for each replicate.", type=int, default=-1);	
	parser.add_argument("-r", dest="num_reps", help="The number of times to perform the randomization test.", type=int, default=1);
	parser.add_argument("-d", dest="conv_divergence", help="Check for convergent divergent sites by setting this to 1.", type=int, default=0);
	parser.add_argument("-t", dest="num_threads", help="Multiple random replicates can be run in parallel. This sets the number of threads to be used. The quotient of -c and -p MUST be a whole number. Default: 1", type=int, default=1);
	parser.add_argument("-o", dest="output_file", help="Output file name for convergent genes/sites.");

	args = parser.parse_args();

	if errorflag == 0:
		if args.input_directory == None or args.output_file == None or args.num_spec == None:
			core.errorOut(1, "Input (-i), output (-o), and the number of species (-s) must all be specified");
			optParse(1);

		if args.num_reps <= 0:
			core.errorOut(2, "The number of replicates (-r) must be a positive integer");
			optParse(1);

		if args.conv_divergence not in [0,1]:
			core.errorOut(3, "-d must take values of either 0 or 1");
			optParse(1);

		if args.num_threads <= 0:
			core.errorOut(4, "The number of threads (-t) must be a positive integer");
			optParse(1);

		if (args.num_reps % args.num_threads) != 0:
				core.errorOut(5, "The quotient of -c and -p MUST be a whole number");
				optParse(1);

		return args.input_directory, args.num_spec, args.num_reps, args.conv_divergence, args.num_threads, args.output_file;

	elif errorflag == 1:
		parser.print_help();
		print()
		sys.exit();

#######

def remGapMiss(oldlist):
	newlist = [];
	for aa in oldlist:
		if aa in ["-", "X"]:
			continue;
		newlist.append(aa);
	return newlist;

#######
def convCheck(cur_c, c, number_specs, d, ins, outs):
#	cur_c = 0;
	init_c = cur_c+1;
	while cur_c < c:
		#if c > 1:

		spec_list = all_specs.values();
		rep_specs = [];
		while len(rep_specs) < number_specs:
			r = random.choice(spec_list);
			rep_specs.append(r);
			spec_list.remove(r);

		outfilename = outs + "_" + str(cur_c+1) + ".txt";
		outfile = open(outfilename, "w");
		outfile.write("# ==============================================================================================\n");
		outfile.write("# \t\t\tConvergence testing\n");
		outfile.write("# \t\t\t" + core.getDateTime() + "\n");
		outfile.write("# Using alignments in:\t\t" + indir + "\n");
		outfile.write("# Randomly choosing " + str(number_specs) + " species and performing " + str(c) + " replicate tests for convergence.\n");
		outfile.write("# This is replicate number " + str(cur_c+1) + "\n");
		outfile.write("# Writing output to:\t\t\t" + outfilename + "\n");
		if d == 0:
			outfile.write("# Checking for convergent sites.\n");
		elif d == 1:
			outfile.write("# Checking for divergent sites.\n");
		outfile.write("# Using species:\t" + ",".join(rep_specs));
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
		outfile.write("# Site#\tTargetSpecs\tChromosome\tGeneID\tAlignLen\tPosition\tTargetAlleles\tBackgroundAlleles\n");
		for align in aligns:
			#numbars, donepercent = core.loadingBar(i, count, donepercent, numbars);
			i = i + 1;
			if align.find(".fa") == -1:
				continue;

			infilename = os.path.join(ins, align);
			gid = "_".join(align.split("_")[:2]);
			chrome = align[align.find("chr"):align.find("chr")+4]

			inseqs = core.fastaGetDict(infilename);

			for t1 in rep_specs:
				for t2 in rep_specs:
					if t1 == t2:
						continue;

					targets = [t1, t2];
					backgrounds = [spec for spec in rep_specs if spec not in targets];

					num_targets_present = 0;
					num_bg_present = 0;
					for title in inseqs:
						if any(t in title for t in targets):
							num_targets_present = num_targets_present + 1;
						if any(b in title for b in backgrounds):
							num_bg_present = num_bg_present + 1;

					if num_targets_present == len(targets) and num_bg_present == len(backgrounds):
						# print "The following gene has all target and background species and will be checked:\t\t" + gid;
						totgenes = totgenes + 1;

						seqlen = len(inseqs[inseqs.keys()[0]]);
						# print "Alignment length\t\t", seqlen;

						t_alleles = {};
						b_alleles = {};

						for x in xrange(len(inseqs[inseqs.keys()[0]])):
							for title in inseqs:
								cur_spec = title[1:].replace("\n","");
								if cur_spec in targets:
									t_alleles[cur_spec] = inseqs[title][x];
								if cur_spec in backgrounds:
									b_alleles[cur_spec] = inseqs[title][x];

							t_states = t_alleles.values();
							#t_gap = t_states.count("-");
							#t_missing = t_states.count("X");
							#t_stop = t_states.count("*");

							b_states = b_alleles.values();
							#b_gap = b_states.count("-");
							#b_missing = b_states.count("X");
							#b_stop = b_states.count("*");

							t_final = remGapMiss(t_states);
							b_final = remGapMiss(b_states);

							if t_final == [] or b_final == []:
								continue;

							if d == 0:
								if len(t_final) == len(targets) and len(b_final) == len(backgrounds) and t_final.count(t_final[0]) == len(t_final) and t_final[0] not in b_final:
									numsites = numsites + 1;
									print core.getTime() + " Convergent site found!";
									print "Filename:\t\t" + align;
									print "Chromosome:\t\t" + chrome;
									print "Gene ID:\t\t" + gid;
									print "Alignment length\t", seqlen;
									print "Target alleles:\t\t" + "".join(t_final);
									print "Background alleles:\t" + "".join(b_final);
									print "---------------";
									outline = str(numsites) + "\t" + ",".join(targets) + "\t" + chrome + "\t" + gid + "\t" + str(seqlen) + "\t" + str(x+1) + "\t" + "".join(t_final) + "\t" + "".join(b_final) + "\n";
									outfile.write(outline);

							elif d == 1:
								if len(t_final) == len(targets) and len(b_final) == len(backgrounds) and t_final.count(t_final[0]) != len(t_final) and b_final.count(b_final[0]) == len(b_final):
									if not any(t in b_final for t in t_final):
										numsites = numsites + 1;
										# print "\nDivergent site found!";
										# print "Filename:\t\t" + align;
										# print "Chromosome:\t\t" + chrome;
										# print "Gene ID:\t\t" + gid;
										# print "Alignment length\t", seqlen;
										# print t_final;
										# print b_final;
										outline = str(numsites) + "\t" + ",".join(targets) + "\t" + chrome + "\t" + gid + "\t" + str(seqlen) + "\t" + str(x+1) + "\t" + "".join(t_final) + "\t" + "".join(b_final) + "\n";
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

indir, specs, reps, cd, num_t, outfix = optParse(0);
all_specs = {"Human":"hg19","Chimp":"panTro4","Gorilla":"gorGor3","Orangutan":"ponAbe2","Gibbon":"nomLeu3","Rhesus":"rheMac3","Crab-eating macaque":"macFas5","Baboon":"papHam1","Green monkey":"chlSab1","Marmoset":"calJac3","Squirrel monkey":"saiBol1","Bushbaby":"otoGar3","Chinese tree shrew":"tupChi1","Squirrel":"speTri2","Lesser Egyptian jerboa":"jacJac1","Prairie hamster":"micOch1","Chinese hamster":"criGri1","Golden hamster":"mesAur1","Mouse":"mm10","Rat":"rn5","Naked mole-rat":"hetGla2","Guinea pig":"cavPor3","Chinchilla":"chiLan1","Brush-tailed rat":"octDeg1","Rabbit":"oryCun2","Pika":"ochPri3","Pig":"susScr3","Alpaca":"vicPac2","Bactrian camel":"camFer1","Dolphin":"turTru2","Killer whale":"orcOrc1","Tibetan antelope":"panHod1","Cow":"bosTau7","Sheep":"oviAri3","Domestic goat":"capHir1","Horse":"equCab2","White rhinocerous":"cerSim1","Cat":"felCat5","Dog":"canFam3","Ferret":"musFur1","Panda":"ailMel1","Pacific walrus":"odoRosDiv1","Weddell seal":"lepWed1","Black flying-fox":"pteAle1","Megabat":"pteVam1","David's myotis bat":"myoDav1","Microbat":"myoLuc2","Big brown bat":"eptFus1","Hedgehog":"eriEur2","Shrew":"sorAra2","Star-nosed mole":"conCri1","Elephant":"loxAfr3","Cape elephant shrew":"eleEdw1","Manatee":"triMan1","Cape golden mole":"chrAsi1","Tenrec":"echTel2","Aardvark":"oryAfe1","Opossum":"monDom5","Platypus":"ornAna1"};

print "==============================================================================================";
print "\t\tPairwise Unique Convergence testing";
print "\t\t\t" + core.getDateTime();
print "Using alignments in:\t\t" + indir
print "Choosing " + str(specs) + " species randomly and performing " + str(reps) + " replicate tests.";
print "The task will be split into " + str(num_t) + " processes and species will be chosen as the processes are split.";
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
		print core.getTime() + " Call " + str(x+1) + ":  convCheck(" + str(start_rep) + " , " + str(stop_rep) + " , " + str(specs) + " , " + str(cd) + ", " + indir + " , " + outfix + ")";
		processes.append(mp.Process(target=convCheck,args=(start_rep,stop_rep,specs,cd,indir,outfix)));
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
	convCheck(0,1,specs,cd,indir,outfix);
print core.getTime() + " Done!";
print "==============================================================================================";
