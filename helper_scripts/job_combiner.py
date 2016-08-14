import sys,os
##################################
specnum = sys.argv[1];
repnum = sys.argv[2];

indir = "/N/dc2/scratch/grthomas/marine_june/anc_reps_proper/" + specnum + "spec/" + repnum + ".job/";
outdir = indir + "codeml_combined/";
codeml_ancdir = outdir + "anc_seqs_fa/";
codeml_outdir = outdir + "codeml_out/";

if indir[len(indir)-1] != "/":
	indir = indir + "/";
if outdir[len(outdir)-1] != "/":
	outdir = outdir + "/";

if not os.path.exists(outdir):
	print "+Creating output directory";
	os.system("mkdir " + outdir);
if not os.path.exists(codeml_outdir):
	print "+Creating codeml_out directory";
	os.system("mkdir " + codeml_outdir);
if not os.path.exists(codeml_ancdir):
	print "+Creating codeml_anc_fa directory";
	os.system("mkdir " + codeml_ancdir);

dirlist = os.listdir(indir);

print "-----------------------------------------";
for jobdir in dirlist:
	if jobdir == ".DS_Store" or jobdir.find("codeml_combined") != -1 or jobdir == "job_submit.sh":
		continue;
#	jobnum = jobdir[jobdir.index(".")+1:];
	fulldir = indir + jobdir + "/";
	print "Current job directory:" + fulldir;
	fulllist = os.listdir(fulldir);

	for adir in fulllist:
		if adir.find("run_codeml") == -1:
			continue;

		codemldir = fulldir + adir + "/";

		print "Current codeml directory: " + codemldir;

		cur_outdir = codemldir + "codeml_out/";
		cur_ancdir = codemldir + "anc_seqs_fa/";

		print "Copying codeml output files...";
		cp_cmd = "cp -r " + cur_outdir + "* " + codeml_outdir;
		print cp_cmd;
		os.system(cp_cmd);

		print "Copying ancestral sequences...";
		cp_cmd = cp_cmd = "cp " + cur_ancdir + "* " + codeml_ancdir;
		print cp_cmd;
		os.system(cp_cmd);


	print "-----------------------------------------";
print "Done!";
