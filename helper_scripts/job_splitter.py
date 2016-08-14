import sys, os, argparse
import core

############################################
#Function Definitions
############################################
def optParse(errflag):
#This function handles the command line options.

	parser = argparse.ArgumentParser(description="Splits a directory of files for running as multiple jobs on a HPC cluster.");

	parser.add_argument("-i", dest="input", help="Input directory containing all files on which your job will be run.");
	parser.add_argument("-s", dest="split_size", help="The # of files to place in each split. Default: 400", type=int, default=400);
	parser.add_argument("-t", dest="job_time", help="The job walltime in hours. Should be entered as simply an integer. Default: 60", type=int, default=60);
	parser.add_argument("-n", dest="job_name", help="The prefix for the job name.");
	parser.add_argument("-o", dest="output", help="An output directory in which the split directories will be placed. Default: [input directory].job/", default="");

	args = parser.parse_args();

	if errflag == 0:

		if args.input == None or args.job_name == None:
			parser.print_help();
			sys.exit();

		if args.split_size <= 0:
			core.errorOut(1, "-s can only take positive, non-zero values");
			optParse(1);

		if args.job_time <= 0:
			core.errorOut(2, "-t can only take positive, non-zero values");
			optParse(1);
		
		return args.input, args.split_size, args.job_time, args.job_name, args.output;

	else:
		parser.print_help();
		print "Exiting program.";
		sys.exit();

##########
def makeScript(joboutdir,x,jt,jp):

		print "Job # " + str(jobnum) + ": Writing job script...";
		jobname = jp + "." + str(x)
		scriptname = joboutdir + jobname + ".pbs";
		outfile = open(scriptname, "w");
		outfile.write("#!/bin/bash\n");
		outfile.write("#PBS -k o\n");
		outfile.write("#PBS -l nodes=1:ppn=1,walltime=" + str(jt) + ":00:00,vmem=8GB\n");
		outfile.write("#PBS -M grthomas@indiana.edu\n");
		outfile.write("#PBS -m abe\n");
		outfile.write("#PBS -N " + jobname + "\n");
		outfile.write("#PBS -j oe\n");
		outfile.write("#PBS -o " + joboutdir + "\n");
		outfile.write("#PBS -d " + joboutdir + "\n");
		outfile.write("module load python/2.7.9\n");
		outfile.write("module load python-modules\n");

		cmd = "time -p python2.7 /N/u/grthomas/Mason/bin/gwct/gwct_codeml.py -i " + joboutdir + " -t /N/dc2/scratch/grthomas/marine_june/multiz_mm59.tre -s 2 -p 1 -a 1 -c /N/u/grthomas/Karst/bin/grand-conv/";

		outfile.write(cmd);

		outfile.close();

		return scriptname;

############################################
#Main Block
############################################

indir, s, jobtime, jobprefix, outdir = optParse(0);

if not os.path.isdir(indir):
	core.errorOut(3, "-i must be a valid directory path");
	optParse(1);
else:
	indir = os.path.abspath(indir) + "/";

if outdir == "":
	outdir = indir[:len(indir)-1] + ".job/";
elif not os.path.isdir(outdir):
	core.errorOut(4, "-o must be a valid directory path");
	optParse(1);

print "=======================================================================";
print "INPUT  | Splitting files in:\t\t" + indir;
print "INFO   | Number of files per split:\t" + str(s);
print "OUTPUT | Output directory:\t\t" + outdir;
print "-------------------------------------";

if not os.path.exists(outdir):
	print "+Creating output directory.";
	os.system("mkdir " + outdir);

filelist = os.listdir(indir);
if ".DS_Store" in filelist:
	filelist.remove(".DS_Store");

numfiles = 0;
numjobs = (len(filelist) / s) + 1;
jobfiles = [];

i = 1;
jobnum = 1;
print "==========";
print "Job # " + str(jobnum) + ": Copying files...";
for each in filelist:		
	if each.find(".fa") == -1:
		continue;
	if i > s:
		jobfiles.append(makeScript(jobdir,jobnum,jobtime,jobprefix));
		jobnum = jobnum + 1;
		print "==========";
		print "Job # " + str(jobnum) + ": Copying files...";
		i = 1;

	jobdir = outdir + str(jobnum) + "/";
	if not os.path.exists(jobdir):
		#print "Job # " + str(jobnum);
		print "+Creating job directory.";
		os.system("mkdir " + jobdir);

	inname = indir + each;
	destname = jobdir + each;
	cp_cmd = "cp " + inname +  " " + destname;
	os.system(cp_cmd);

	i = i + 1;

jobfiles.append(makeScript(jobdir,jobnum,jobtime,jobprefix));
print "==========";
print "Writing job_submit script...";
bashscript = outdir + "job_submit.sh";
bfile = open(bashscript, "w");
bfile.write("#!/bin/bash\n");
for each in jobfiles:
	bfile.write("qsub " + each + "\n");
bfile.close();
print "==========";
print "Done!";
print "=======================================================================";
