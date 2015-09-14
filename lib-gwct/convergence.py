#############################################################################
#This is the main function of the GWCT. It takes a list of targets, a gene,
#and a tree as the main inputs and detects convergent, divergent, and unique
#substitutions from a sequence file with ancestral states.
#Gregg Thomas
#September 2015.
#############################################################################

import sys
import gwctree
import gwctlib

#############################################################################
def convCheck(in_name, c_name, d_name, target_list, c_sites, d_sites, p_thresh, chrome, g, td):

#	print treefilename
#	print infilename
#	print tree;
#	print td;
#	print newtree;
	#Prep info for this gene and reading of the tree file.

#	print "-----";
#	print gid;

	for t in xrange(len(target_list)):
		if type(target_list[t]) == str:
			if "_" in target_list[t]:
				target_list[t] = gwctree.specRelabel(target_list[t],td);
		else:
			for s in xrange(len(target_list[t])):
				if "_" in s:
					target_list[t][s] = gwctree.specRelabel(target_list[t][s],td);
	#Relabeling of target labels to match those in the tree.

	conv_nodes = {};
	for t in target_list:
		if type(t) == list:
			f, com_anc = gwctree.comAnc(t,td);
			conv_nodes[">node #" + com_anc] = ">node #" + td[com_anc][1];
		else:
			if "_" in t:
				conv_nodes[">" + t[t.index("_")+1:]] = ">node #" + td[t][1];
			else:
				conv_nodes[">node #" + t] = ">node #" + td[t][1];
	#Identification of ancestral and target nodes in the tree. If more than one species is present in a group, this will get the common ancestor
	#of those species as the target node.

	probseqs = gwctlib.fastaGetDict(in_name);
	#Reading of the sequence probability file.

	probs = {};
	seqs = {};
	for seq in probseqs:
		if seq.find(" ") != -1:
			probseqs[seq] = filter(None, probseqs[seq].split(" "));
		#print len(probseqs[seq]);
	seqlen = len(probseqs.values()[0]);

	for seq in probseqs:
		if type(probseqs[seq]) == str:
			seqs[seq] = probseqs[seq];
			probs[seq] = ['1.0'] * seqlen;
		else:
			seqs[seq] = "";
			probs[seq] = [];
			for pos in probseqs[seq]:
				bp = pos.split("_");
				b = bp[0];
				p = bp[1];
				seqs[seq] = seqs[seq] + b;
				probs[seq].append(p);
	#Splitting of the sequences and probabilities into separate dictionaries.


	x = 0;
	while x < seqlen:
	#For each site in the current gene, check for convergence. x is the index of the current residue.
		#site = {};
		#for seq in seqs:
		#	site[seq] = seqs[seq][x];
		anc_site = {};
		conv_site = {};
		cur_probs = [];
		for seq in seqs:
			if seq in conv_nodes.keys():
				conv_site[seq] = seqs[seq][x];
				cur_probs.append(float(probs[seq][x]));
			elif seq in conv_nodes.values():
				cur_probs.append(float(probs[seq][x]));
				anc_site[seq] = seqs[seq][x];
		#Retrieves states and probabilities for target and ancestral nodes.
		#print conv_nodes;
		#print conv_site;
		#print anc_site;
		#print "------";
		#sys.exit();

		if "X" in conv_site.values() or "-" in conv_site.values():
			x = x + 1;
			continue;

		#if p_thresh > 0 and all(c >= p_thresh for c in cur_probs):
		if all(c >= p_thresh for c in cur_probs):
			if conv_site.values().count(conv_site.values()[0]) == len(conv_site.values()) and conv_site.values()[0] not in anc_site.values():
				outline = str(c_sites+1) + "\t" + chrome + "\t" + g + "\t" + str(seqlen) + "\t" + str(x+1) + "\t" + "".join(anc_site.values()) + "\t" + "".join(conv_site.values()) + "\n";
				convfile = open(c_name, "a");

				c_sites = c_sites + 1;
				convfile.write(outline);
				convfile.close();

				#print "Convergent site found!!";
				#print c_name;
				#print target_list;			
				#print conv_site;
				#print anc_site;
				#print cur_probs;
			#If a convergent site is found, write it to the file and increment c_sites. Also checkes if a probability threshold has been set and only writes if all probs for all
			#reconstructed states are >= to that threshold.

			if all(anc_site[conv_nodes[node]] != conv_site[node] for node in conv_nodes) and all(conv_site.values().count(aa) == 1 for aa in conv_site.values()):
				outline = str(c_sites+1) + "\t" + chrome + "\t" + g + "\t" + str(seqlen) + "\t" + str(x+1) + "\t" + "".join(anc_site.values()) + "\t" + "".join(conv_site.values()) + "\n";
				divfile = open(d_name, "a");

				d_sites = d_sites + 1;
				divfile.write(outline);
				divfile.close();

				#print "Divergent site found!!";
				#print d_name;
				#print target_list;
				#print conv_nodes;		
				#print conv_site;
				#print anc_site;
				#print cur_probs;
				#sys.exit();

		x = x + 1;

	return c_sites, d_sites;

#############################################################################
