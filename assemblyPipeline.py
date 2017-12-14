'''This script will call FastQC, Trimmomatic (if necessary), ABYSS, 
and filterMetagenomicSequences.py.'''

import os
from subprocess import Popen
from shlex import split
from glob import glob
from zipfile import ZipFile
import shutil

DEVNULL = open(os.devnull, "w")

def minLen(inset, k):
	# Set minimum length to K if it is not given
	" SLIDINGWINDOW:4:25 LEADING:28 TRAILING:28 MINLEN:"
	i = inset.find("MINLEN:")
	mlen = inset[i:]
	if " " in mlen:
		mlen = mlen[:mlen.find(" ")]
	if mlen.split(":")[1]:
		# Return orignal settings if minlen is specified
		return inset
	else:
		# Set minlen to k
		return inset + k 

def fastQC(fastqs, outpath, cpu, trim=False):
	# Calls FastQC for all input files
	print("\n\tRunning FastQC...\n")
	quals = []
	cmd = "perl /opt/FastQC/fastqc -t " + cpu + " -o "
	for batch in fastqs:
		# Construct cammand and output directory for each batch
		outdir = outpath + "fastqc-" + batch + "/"
		if not os.path.isdir(outdir):
			os.mkdir(outdir)
		string = cmd +  outdir
		for sample in fastqs[batch]:
			# Add paths to each file
			string += (string + " " + fastqs[batch][sample][0] + 
						" " + fastqs[batch][sample][1])
		# Call FastQC, save output directory
		if trim == False:
			try:
				fqc = Popen(split(string))
				fqc.communicate()
				quals.append(outdir)
			except:
				print("\tFastQC returned an error.")
				quit()
		else:
			quals.append(outdir)
	return quals

def readSum(archive, summ):
	# Checks summary file for adapter content and sequence quality
	# Ignores warnings by default 
	squal = True
	adapter = False
	with archive.open(summ) as s:
		summary = s.readlines()
	if ("FAIL").encode() in summary[1]:
		squal = False
	if ("FAIL").encode() in summary[10]:
		adapter = True
	return adapter, squal

def parseQuality(quals, outpath):
	# Checks FastQC output and determines if trimming is necessary
	print("\n\tExamining FastQC reports...")
	afails = []
	trim = False
	for i in quals:
		reports = glob(i + "*_fastqc.zip")
		for rep in reports:
			# Isolate archive name 
			summ = rep[rep.rfind("/")+1:].replace(".zip", "") + "/summary.txt"
			archive = ZipFile(rep)
			adapter, squal = readSum(archive, summ)
			if adapter == True:
				# Add file with addpaters to fail list
				afails.append(rep[rep.rfind("/")+1:rep.rfind("_")])
			elif squal == False:
				# Trim all files if any fail quality check
				trim = True					
	if afails:
		print("\tTrim adapters from the following files:")
		for i in afails:
			print("\t\t" + i)
	return trim

#-----------------------------------------------------------------------------

def trimmomatic(fastqs, outpath, cpu, trm, notrim=False):
	# Calls Trimmomatic on fastq files and edits fastq paths
	print("\n\tRunning Trimmomatic...\n")
	paired = {}
	flist = []
	cmd1 = ("java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -threads " 
			+ cpu + " -trimlog ")
	for batch in fastqs:
		# Construct cammand and output directory for each batch
		paired[batch] = {}
		unpaired = outpath + "unpaired-" + batch + "/"
		paired = outpath + "paired-" + batch + "/"
		log = paired + "trimLog.txt"
		for i in [paired, unpaired]:
			if not os.path.isdir(i):
				os.mkdir(i)
		for sample in fastqs[batch]:
			# Construct command for each pair of reads
			r1 = fastqs[batch][sample][0]
			r2 = fastqs[batch][sample][1]
			# Isolate file names
			n1 = r1[r1.rfind("/")+1:]
			n2 = r2[r2.rfind("/")+1:]
			# Add new paths
			paired[batch][sample] = [paired + n1, paired + n2]
			string = (log + " " + r1 + " " + r2 + " " + paired + n1 + " " 
					+ unpaired + n1 + " " + paired + n2 + " " + unpaired + n2)
			if notrim == False:
				# Call Trimmomatic
				try:
					tmm = Popen(split(cmd1 + string +  " " + trm))
					tmm.communicate()
					# Remove unpaired directory
					shutil.rmtree(unpaired)
				except:
					flist.append(sample + " " + n1)
					flist.append(sample + " " + n2)
	if flist:
		print("\nThe following files failed trimming:")
		for i in flist:
			print("\t" + i)
		quit()
	return paired, True

def abyss(fastqs, outpath, cpu, k, startdir):
	# Calls ABYSS for each batch
	print("\n\tRunning ABYSS assembler...")
	contigs = {}
	for batch in fastqs:
		# Construct cammand and output directory for each batch
		outdir = outpath + "abyss-" + batch + "/"
		if not os.path.isdir(outdir):
			os.mkdir(outdir)
		os.chdir(outdir)
		cmd = "/usr/bin/abyss-pe np=" + cpu + " k=" + k + " name=" + batch
		if len(fastqs[batch]) == 1:
			# Construct cammand for single library
			cmd += " in='"
			for sample in fastqs[batch]:
				s = fastqs[batch][sample]
				cmd += s[0] + " " + s[1] + "'"
		else:
			# Construct command for multiple libraries
			lib = ""
			libs = []
			for sample in fastqs[batch]:
				lib += sample + " "
				s = fastqs[batch][sample]
				libs.append(sample + "='" + s[0] + " " + s[1] + "'")
			# Append library and file pairs to command
			cmd += " lib='" + lib[:-1] + "'"
			for i in libs:
				cmd += " " + i
			cmd += "'"
		try:
			ab = Popen(split(cmd))
			ab.communicate()
			# Append path to output file
			contigs[batch] = outdir + batch + "-8.fa"
		except:
			print(("\tError: {} failed ABYSS assembly.").format(batch))
	os.chdir(startdir)
	return contigs

def metaSPAdes(fastqs, outpath, cpu):
	# Calls MetaSPAdes fro each pair of reads
	print("\n\tRunning SPAdes in metagomic assembly mode...")
	contigs = {}
	for batch in fastqs:
		if len(fastqs[batch]) > 1:
			print("\tError. MetaSPAdes can only run on single paired-end \
libraries. Make sure there is only one pair of samples per batch.\n\tExiting.\n")
			quit()
		# Construct cammand and output directory for each batch
		outdir = outpath + "spades-" + batch + "/"
		if not os.path.isdir(outdir):
			os.mkdir(outdir)
		for sample in fastqs[batch]:
			# Get file paths
			s = fastqs[batch][sample]
		cmd = ("python /opt/SPAdes/bin/spades.py --meta -t {} --pe1-1 {} --pe1-2 {} \
-o {}").format(cpu, s[0], s[1], outdir)
		try:
			sp = Popen(split(cmd))
			sp.comunicate()
			# Append path to output file
			contigs[batch] = outdir + "contigs.fasta"
		except:
			print(("\tError. {} failed SPAdes assembly.").format(batch))
	return contigs
