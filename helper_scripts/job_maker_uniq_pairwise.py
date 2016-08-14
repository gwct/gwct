#11 - 57

x = 15;
while x <= 55:
	if x % 5 == 0:
		print x;
		outfilename = "/N/dc2/scratch/grthomas/marine_june/uniq_reps_pairwise/job_scripts/uniq_pw_" + str(x) + ".pbs";
		outfile = open(outfilename, "w");
		outfile.write("#!/bin/bash\n");
		outfile.write("#PBS -k o\n");
		outfile.write("#PBS -l nodes=1:ppn=6,walltime=100:00:00,vmem=12GB\n");
		outfile.write("#PBS -M grthomas@indiana.edu\n");
		outfile.write("#PBS -m abe\n");
		outfile.write("#PBS -N mm.uniq.pw."+str(x)+"\n");
		outfile.write("#PBS -j oe\n");
		outfile.write("#PBS -o /N/dc2/scratch/grthomas/marine_june/uniq_reps_pairwise/\n");
		outfile.write("#PBS -d /N/dc2/scratch/grthomas/marine_june/uniq_reps_pairwise/\n");
		# outfile.write("module load python/2.7.9\n");
		# outfile.write("module load python-modules\n");

		cmd = "time -p python2.7 /N/dc2/scratch/grthomas/marine_june/scripts/convergence_checker_mp_pairwise.py -i /N/dc2/scratch/grthomas/marine_june/59spec_src/ -s " + str(x) + " -r 5 -o /N/dc2/scratch/grthomas/marine_june/scripts/test/ -t 5";
		outfile.write(cmd);

		outfile.close();

	x = x + 1;

