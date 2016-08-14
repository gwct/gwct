import sys, os, random
import core

numspec = int(sys.argv[1]);
rep = sys.argv[2];


#indir = "/Users/Gregg/Desktop/marine_2/seq_files/conv_reps_full/59spec/";
#outdir = "/Users/Gregg/Desktop/marine_2/seq_files/june/anc_reps_proper/"+str(numspec)+"spec/"+rep+"/";

indir = "/N/dc2/scratch/grthomas/marine_june/anc_reps_proper/59spec_src/";
outdir = "/N/dc2/scratch/grthomas/marine_june/anc_reps_proper/"+str(numspec)+"spec/"+rep+"/";

#print outdir;

if not os.path.exists(outdir):
	print "+Creating output directory";
	print outdir;
	print "------";
	os.system("mkdir " + outdir);

print "+Selecting random background species"
spec_to_get = [">triMan1", ">odoRosDiv1", ">lepWed1", ">orcOrc1", ">turTru2"];
possible_spec = [">hg19",">panTro4",">gorGor3",">ponAbe2",">nomLeu3",">rheMac3",">macFas5",">papHam1",">chlSab1",">calJac3",">saiBol1",">otoGar3",">tupChi1",">speTri2",">jacJac1",">micOch1",">criGri1",">mesAur1",">mm10",">rn5",">hetGla2",">cavPor3",">chiLan1",">octDeg1",">oryCun2",">ochPri3",">susScr3",">vicPac2",">camFer1",">panHod1",">bosTau7",">oviAri3",">capHir1",">equCab2",">cerSim1",">felCat5",">canFam3",">musFur1",">ailMel1",">pteAle1",">pteVam1",">myoDav1",">myoLuc2",">eptFus1",">eriEur2",">sorAra2",">conCri1",">loxAfr3",">eleEdw1",">chrAsi1",">echTel2",">oryAfe1",">monDom5",">ornAna1"]
while len(spec_to_get) < numspec:
	r = random.choice(possible_spec);
	spec_to_get.append(r);
	possible_spec.remove(r);
print "------";
print "+Getting sequences"
filelist = os.listdir(indir);
numgenes = 0;

for each in filelist:
	if each.find(".fa") == -1:
		continue;

	infilename = indir + each;
	inseqs = core.fastaGetDict(infilename);

	if all(s in inseqs for s in spec_to_get):
		#print "Getting: " + each;
		numgenes = numgenes + 1;
		outfilename = outdir + each;
		ofile = open(outfilename, "w");
		ofile.write("");
		ofile.close();
		for seq in inseqs:
			if seq in spec_to_get:
				core.writeSeqOL(outfilename, inseqs[seq], seq);

print "Got " + str(numgenes) + " genes.";
print "-----------";

