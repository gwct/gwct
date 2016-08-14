#11 - 57

x = 11;
while x < 55:
	print x;
	outfilename = "/N/dc2/scratch/grthomas/marine_june/uniq_reps_proper/job_scripts/uniq_" + str(x) + ".pbs";
	outfile = open(outfilename, "w");
	outfile.write("#!/bin/bash\n");
	outfile.write("#PBS -k o\n");
	outfile.write("#PBS -l nodes=1:ppn=10,walltime=8:00:00,vmem=8GB\n");
	outfile.write("#PBS -M grthomas@indiana.edu\n");
	outfile.write("#PBS -m abe\n");
	outfile.write("#PBS -N mm.uniq."+str(x)+"\n");
	outfile.write("#PBS -j oe\n");
	outfile.write("#PBS -o /N/dc2/scratch/grthomas/marine_june/uniq_reps_proper/\n");
	outfile.write("#PBS -d /N/dc2/scratch/grthomas/marine_june/uniq_reps_proper/\n");
	# outfile.write("module load python/2.7.9\n");
	# outfile.write("module load python-modules\n");

	cmd = "time -p python2.7 /N/dc2/scratch/grthomas/marine_june/scripts/convergence_checker_mp.py -i /N/dc2/scratch/grthomas/marine_june/59spec_src/ -t 'odoRosDiv1 lepWed1 turTru2 orcOrc1 triMan1' -r "+str(x)+" -c 100 -p 10 -o /N/dc2/scratch/grthomas/marine_june/uniq_reps_proper/"+str(x)+"spec_100reps/"+str(x)+"spec"
	outfile.write(cmd);

	outfile.close();

	x = x + 1;
