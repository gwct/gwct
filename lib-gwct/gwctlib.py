#############################################################################
#A library of functions to help out with the GWCT.
#Functions taken from CORE repository (www.github.com/gwct/core.git)
#Gregg Thomas
#Forked from CORE September 2015
#############################################################################

import string
import sys
import os
import re
import subprocess
import datetime

#############################################################################

def loadingBar(counter, length, done, bars):
#This function serves as a text loading bar for long scripts with counters. The following
#lines must be added within the script to initialize and terminate the script:
#Initilization:
#numlines = core.getFileLen(alnfilename);
#numbars = 0;
#donepercent = [];
#i = 0;
#Termination:
#	pstring = "100.0% complete.";
#	sys.stderr.write('\b' * len(pstring) + pstring);
#	print "\nDone!";
#
#If length is lines in a file use the core.getFileLen function to get the number of lines in the file

	percent = float(counter) / float(length) * 100.0;
	percentdone = int(percent);

	p = str(percent)
	pstring = " " + p[:5] + "% complete.";

	if percentdone % 2 == 0 and done != None and percentdone not in done:
		loading = "";
		loading = "[";
		j = 0;
		while j <= bars:
			loading = loading + "*";
			j = j + 1;
		while j < 50:
			loading = loading + "-";
			j = j + 1;
		loading = loading + "]";

		loading = loading + "                 ";
		sys.stderr.write('\b' * len(loading) + loading);

		done.append(percentdone);
		bars = bars + 1;

	sys.stderr.write('\b' * len(pstring) + pstring);

	return bars, done;

#############################################################################

def getDateTime():
#Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y | %I:%M:%S");

#############################################################################

def getTime():
#Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%I:%M:%S");

#############################################################################

def getLogTime():
#Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y-%I.%M.%S");

#############################################################################

def printWrite(o_name, o_line):
#Function to print a string AND write it to the file.
	print o_line;
	f = open(o_name, "a");
	f.write(o_line + "\n");
	f.close();

#############################################################################

def logCheck(lopt, lfilename, outline):
#Function checks whether or not to write to a logfile, print something, or both.
	if lopt == 1:
		printWrite(lfilename, outline);
	else:
		print outline;

#############################################################################

def errorOut(errnum, errmsg):
#Formatting for error messages.
	fullmsg = "|**Error " + str(errnum) + ": " + errmsg + " |";
	border = " " + "-" * (len(fullmsg)-2);
	print border + "\n" + fullmsg + "\n" + border;

#############################################################################

def getOutdir(indir, prefix, stime):
#Retrieves full input directory name and proper output directory name for other scripts.
	if not os.path.isdir(indir):
		errorOut(0, "-i must be a valid directory path");
		sys.exit();
	indir = os.path.abspath(indir);
	if indir[-1] != "/":
		indir = indir + "/";
	filelist = os.listdir(indir);
	used = [];
	for each in filelist:
		if each.find("-" + prefix) != -1:
			used.append(int(each[:each.index("-")]));
	if used != []:
		outdir = indir + str(max(used)+1) + "-" + prefix + "-" + stime + "/";
	else:
		outdir = indir + "1-" + prefix + "-" + stime + "/";

	return indir, outdir;

#############################################################################

def spacedOut(string, totlen):
#Properly adds spaces to the end of a message to make it a given length
	spaces = " " * (totlen - len(string));
	return string + spaces;

#############################################################################

def filePrep(filename, header):
#Writes over a file, header optional (if no header just pass "")
	f = open(filename, "w");
	f.write(header);
	f.close();

#############################################################################
##########################################################################################################################################################
#SEQUENCE FORMAT READERS AND WRITERS
##########################################################################################################################################################
#FASTA
#############################################################################

def fastaGetDict(i_name):
#fastaGetDict reads a FASTA file and returns a dictionary containing all sequences in the file with 
#the key:value format as title:sequence.

	seqdict = {};
	for line in open(i_name, "r"):
		line = line.replace("\n", "");
		if line[:1] == '>':
			curkey = line;
			seqdict[curkey] = "";
		else:
			seqdict[curkey] = seqdict[curkey] + line;

	return seqdict;

#############################################################################

def fastaGetFileInd(i_name):
#fastaGetFileInd reads a FASTA file and returns a dictionary containing file indexes for each title
#and sequence with the key:value format as [title start index]:[sequence start index]

	infile = open(i_name, "r");
	indList = [];
	firstflag = 0;
	curlist = [];

	line = "derp";

	while line != '':
		line = infile.readline();
		if line[:1] == '>':
			if firstflag == 1:
				curseqend = infile.tell() - len(line) - 1;
				curlist.append(curseqend);
				indList.append(curlist);
				curlist = [];

			curtitlestart = infile.tell() - len(line);
			curtitleend = infile.tell() - 1;
			curseqstart = infile.tell();

			curlist.append(curtitlestart);
			curlist.append(curtitleend);
			curlist.append(curseqstart);

			firstflag = 1;

	curseqend = infile.tell() - len(line) - 1;
	curlist.append(curseqend);
	indList.append(curlist);

	infile.close();
	return indList;
		
#############################################################################

def getFastafromInd(i_name, titlestart, titleend, seqstart, seqend):
#This takes the file index for a corresponding FASTA title and sequence (as retrieved by
#fastaGetFileInd and returns the actual text of the title and the sequence.

	infile = open(i_name, "r");

	infile.seek(titlestart);
	title = infile.read(titleend - titlestart);

	infile.seek(seqstart);
	seq = infile.read(seqend - seqstart);

	infile.close();

	title = title.replace("\n", "");
	seq = seq.replace("\n", "");

	return title, seq;

#############################################################################
