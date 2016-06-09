import sys, argparse, core, treeparse

############################################
#Function Definitions
############################################

def optParse(errorflag):
#This function handles the command line options.
	parser = argparse.ArgumentParser();
	parser.add_argument("-i", dest="input_file", help="A file containing a FASTA alignment with ancestral sequences from codeml.");
	parser.add_argument("-t", dest="tree_file", help="The corresponding tree file from the codeml ancestral reconstructions.");
	parser.add_argument("-s", dest="site_num", help="The site in the alignment that you wish to map to the tree.", type=int);
	parser.add_argument("-o", dest="output_file", help="The name of an output file to write the new tree.");
	args = parser.parse_args();

	if errorflag == 0:
		if args.input_file == None or args.tree_file == None or args.site_num == None:
			core.errorOut(1, "-i, -t, and -s must always be defined");
			optParse(1);

		if args.output_file == None:
			of = "";
		else:
			of = args.output_file;

		return args.input_file, args.tree_file, args.site_num, of;

	elif errorflag == 1:
		parser.print_help();
		sys.exit();

############################################
#Main Block
############################################

infilename, treefilename, s, outfilename = optParse(0);

print "=======================================================================";
print "Mapping column", s, "from the alignment file", infilename, "to the tree in", treefilename;
print "-------------------------------------"

tree = open(treefilename, "r").read().replace("\n","");
td, newtree = treeparse.treeParse(tree,2);

inseqs = core.fastaGetDict(infilename);
site = {};
for seq in inseqs:
	if "#" in seq:
		newseq = seq[seq.index("#")+1:];
	else:
		newseq = seq[1:];
	site[newseq] = inseqs[seq][s];

for node in td:
	if "_" in node:
		new_node = node + "_" + site[node[node.index("_")+1:]];
	else:
		new_node = node + "_" + site[node];
	newtree = newtree.replace(node,new_node);
newtree = newtree + ";"
print newtree;
if outfilename != "":
	outfile = open(outfilename, "w");
	outfile.write(newtree);
	outfile.close();

print "-------------------------------------"
print "Done!"
print "=======================================================================";
