#############################################################################
#This is the main function of the GWCT. It takes a list of targets, a gene,
#and a tree as the main inputs and detects convergent, divergent, and unique
#substitutions from a sequence file with ancestral states.
#Gregg Thomas
#September 2015.
#############################################################################

import sys, gwctree, gwctcore

#############################################################################
def convCheck(in_name, conv_dict, conv_key, target_list, p_thresh, g, td, pair_opt):

	tmp_list = [];
	target_alleles = {};
	for t in target_list:
		new_targets = [];
		for s in t:
			node = gwctree.specRelabel(s,td);
			if node not in td:
				print "* WARNING: Target node '" + s + "' not in alignment. Removing from target list.";
			else:
				target_alleles[">" + s] = "";
				new_targets.append(node);
		if new_targets != []:
			tmp_list.append(new_targets);
	target_list = tmp_list;
	# Checks to make sure every target node is in the alignment.

	# print 1,conv_key
	# print 2,target_list;
	# print 3,target_alleles;
	# print "-----";

	# for t in xrange(len(target_list)):
	# 	for s in xrange(len(target_list[t])):
	# 		target_alleles[">" + target_list[t][s]] = "";
	# 		target_list[t][s] = gwctree.specRelabel(target_list[t][s],td);
	# Relabeling of target labels to match those in the tree.

	conv_nodes = {};
	for t in target_list:
		if len(t) == 1:
			t = t[0];
			if "_" in t:
				conv_nodes[">" + t[t.index("_")+1:]] = ">node #" + td[td[t][1]][3];
			else:
				conv_nodes[">node #" + td[t][3]] = ">node #" + td[td[t][1]][3];
		else:
			f, com_anc = gwctree.comAnc(t,td);
			conv_nodes[">node #" + td[com_anc][3]] = ">node #" + td[td[com_anc][1]][3];
	# Identification of ancestral and target nodes in the tree. If more than one species is present in a group, this will get the common ancestor
	# of those species as the target node.
	# print conv_nodes;
	# print "-----";

	if pair_opt == 1:
		target_alleles = { ">" + t[t.index("_")+1:] : '' for t in target_alleles };

	probseqs = gwctcore.fastaGetDict(in_name);
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
		bg_alleles = {};

		for seq in seqs:
			if seq in conv_nodes.keys():
				conv_site[seq] = seqs[seq][x];
				cur_probs.append(float(probs[seq][x]));
			elif seq in conv_nodes.values():
				cur_probs.append(float(probs[seq][x]));
				anc_site[seq] = seqs[seq][x];
			if seq.find("node #") == -1 and seq not in target_alleles.keys():
				bg_alleles[seq] = seqs[seq][x];
			elif seq.find("node #") == -1:
				target_alleles[seq] = seqs[seq][x];

		#Retrieves states and probabilities for target and ancestral nodes.
		#print target_alleles;
		#print conv_nodes;
		#print conv_site;
		#print anc_site;
		#print bg_alleles;
		#print len(bg_alleles);
		#print len(seqs.keys());
		#print "------";
		#sys.exit();

		if "X" in conv_site.values() or "-" in conv_site.values():
			x = x + 1;
			continue;

		#if p_thresh > 0 and all(c >= p_thresh for c in cur_probs):
		if all(c >= p_thresh for c in cur_probs):
			if conv_site.values().count(conv_site.values()[0]) == len(conv_site.values()) and conv_site.values()[0] not in anc_site.values():
				#c_sites = c_sites + 1;
				outline = g + "\t" + str(seqlen) + "\t" + str(x+1) + "\t" + "".join(anc_site.values()) + "\t" + "".join(conv_site.values()) + "\n";
				conv_dict[conv_key][0].append(outline);
				#c_sites.append(outline);
				#convfile = open(c_name, "a");

				#convfile.write(outline);
				#convfile.close();

				#print "Convergent site found!!";
				#print c_name;
				#print target_list;
				#print conv_site;
				#print anc_site;
				#print cur_probs;
			#If a convergent site is found, write it to the file and increment c_sites. Also checks if a probability threshold has been set and only writes if all probs for all
			#reconstructed states are >= to that threshold.

			if all(anc_site[conv_nodes[node]] != conv_site[node] for node in conv_nodes) and all(conv_site.values().count(aa) == 1 for aa in conv_site.values()):
				#d_sites = d_sites + 1;
				outline = g + "\t" + str(seqlen) + "\t" + str(x+1) + "\t" + "".join(anc_site.values()) + "\t" + "".join(conv_site.values()) + "\n";
				conv_dict[conv_key][1].append(outline);
				#d_sites.append(outline);
				#divfile = open(d_name, "a");

				#divfile.write(outline);
				#divfile.close();

				#print "Divergent site found!!";
				#print d_name;
				#print target_list;
				#print conv_nodes;
				#print conv_site;
				#print anc_site;
				#print cur_probs;
				#sys.exit();

		if pair_opt in [0,1] and target_alleles.values().count(target_alleles.values()[0]) == len(target_alleles.values()) and target_alleles.values()[0] not in bg_alleles.values() and "-" not in bg_alleles.values() and "X" not in bg_alleles.values():
				#u_sites = u_sites + 1;
				outline = g + "\t" + str(seqlen) + "\t" + str(x+1) + "\t" + "".join(bg_alleles.values()) + "\t" + "".join(target_alleles.values()) + "\n";
				conv_dict[conv_key][2].append(outline);
				#u_sites.append(outline);
				#unifile = open(u_name, "a");

				#unifile.write(outline);
				#unifile.close();

				#print "Unique site found!!";
				#print u_name;
				#print target_alleles;
				#print bg_alleles;
				#sys.exit();

		x += 1;

	return conv_dict;

#############################################################################
