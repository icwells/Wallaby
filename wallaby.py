'''This script will call FastQC, Trimmomatic (if necessary), ABYSS, 
and filterMetagenomicSequences.py, and concatenate linear and circular 
contig files for BLAST.'''

import argparse
import os
import sys
import MySQLdb
from getpass import getpass
from datetime import datetime
from multiprocessing import cpu_count
from assemblyPipeline import *
from filterMetagenomicSequences import *

# Import blast scripts from Kiwi
KIWI = "../Kiwi/bin/"
sys.path.append(KIWI)
from blastSeqs import blastSeqs, ublast
from blastResults import readBlast, getSeqs, printOutput

def Ublast(queries, conf):
	# Calls ublast pipeline and concatenates output
	os.chdir(KIWI)
	db = MySQLdb.connect("localhost", conf["name"], conf["pw"], "ASUviralDB")
	for i in queries:
		done = False
		outdir = i[:i.rfind("-")+1] + "Ublast/"
		if not os.path.isdir(outdir):
			os.mkdir(outdir)
		nout = outdir + "ublastN.outfmt6"
		xout = outdir + "ublastX.outfmt6"
		# Call BLAST pipeline
		done = ublast(i, outdir, conf["db"], conf["cpu"])
		if done == True:
			# Concatenate ouput
			blastn, blastx = readBlast(xout, nout)
			blastn, blastx = getSeqs(i, blastn, blastx)
			printOutput(db, outdir, blastn, blastx)
	db.close()

def blast(queries, conf):
	# Calls blast pipeline and concatenates output
	db = MySQLdb.connect("localhost", conf["name"], conf["pw"], "ASUviralDB")
	for i in queries:
		done = False
		outdir = i[:i.rfind("-")+1] + "Blast/"
		if not os.path.isdir(outdir):
			os.mkdir(outdir)
		nout = outdir + "blastn.outfmt6"
		xout = outdir + "blastx.outfmt6"
		# Call BLAST pipeline
		done = blastSeqs(i, outdir, conf["db"], conf["cpu"])
		if done == True:
			# Concatenate ouput
			blastn, blastx = readBlast(xout, nout)
			blastn, blastx = getSeqs(i, blastn, blastx)
			printOutput(db, outdir, blastn, blastx)
	db.close()

def sort(contigs, path, align):
	# Calls filterMetagenomicSequences.py
	queries = []
	for i in contigs:
		# Make outdirs and file prefixes
		name = i[i.rfind("/")+1:i.rfind("-")]
		outdir = path + name + "-Sorted/"
		outpath = outdir + name
		if align == False:
			# Skip sorting if resuming from Blast
			if not os.path.isdir(outdir):
				os.mkdir(outdir)
			fasta = sortFasta(i)
			circ = circularSeqs(fasta, outpath)
			printContigs(outpath, fasta, circ)
		queries.append(outpath + "_min250bp.fasta")
	return queries	

def inputDict(manifest):
	# Saves information from manifest file
	fastqs = {}
	with open(manifest, "r") as mani:
		for line in mani:
			if  line[0] == "#":
				# skip commented lines
				pass
			elif line.strip():
				# Proceed if line is not blank
				line = line.strip().split("\t")
				if len(line) > 3:
					print("\n\tManifest file is malformed.")
					quit()
				batch = line[2]
				sample = line[1]
				# Create dict of fastq files for each sample and add to fastqs
				if batch in fastqs.keys():
					if sample in fastqs[batch].keys():
						fastqs[batch][sample].append(line[0])
					else:
						fastqs[batch][sample] = [line[0]]
				else:
					fastqs[batch] = {}
					fastqs[batch][sample] = [line[0]]
	for batch in fastqs:
		for sample in fastqs[batch]:
			if len(fastqs[batch][sample]) != 2:
				print("\n\tError: Incorect number of reads found in sample:", sample)
				print("\tPlease supply a forward and reverse read for each sample.")
				quit()
	return fastqs

#-----------------------------------------------------------------------------

def config(noblast):
	# Reads config file and returns a dict of settings
	run = True
	# Make empty dict with required settings
	if noblast == False:
		conf = {"cpu":"", "k":"", "trim":"", "db":"", "name":""}
	elif noblast == True:
		conf = {"cpu":0, "k":"", "trim":""}
	if not os.path.isfile("config.txt"):
		# Quit if no config file is present
		print("\tCannot find config.txt. Exiting.")
		quit()
	with open("config.txt", "r") as config:
		for line in config:
			line = line.strip()
			if line:
				if line[0] != "#":
					# Add settings to dict
					line = line.split("=")
					pair = []
					for i in line:
						pair.append(i.strip())
					if pair[0] == "threads":
						conf["cpu"] = pair[1]
						# Check cpu usage
						if int(conf["cpu"]) > cpu_count():
							conf["cpu"] = str(cpu_count())
					elif pair[0] == "k":
						conf["k"] = pair[1]
					elif pair[0] == "trim":
						# Make sure trimmomatic settings a leading space
						conf["trim"] = " " + pair[1]
					elif "databases" in pair:
						if pair[1][-1] != "/":
							conf["db"] = pair[1] + "/"
						else:
							conf["db"] = pair[1]
					elif pair[0] == "username":
						conf["name"] = pair[1]
	for i in conf.keys():
		# Make sure settings are specified
		if not conf[i]:
			print(("\tConfig file missing {}").format(i))
			run = False
	if run == False:
		print("\tExiting.")
		quit()
	else:
		if noblast == False:
			conf["pw"] = getpass(prompt = "\tEnter MySQL password: ")
			try:
				# Test passwoed
				db = MySQLdb.connect("localhost", conf["name"], conf["pw"], "ASUviralDB")
				db.close()
			except:
				print("\nIncorrect password. Access denied.")
				quit()
		return conf

def main():
	starttime = datetime.now()
	startdir = os.getcwd()
	parser = argparse.ArgumentParser(description = "Wallaby version 0.3 (9/14/17) \
script will run FastQC, Trimmomatic (if necessary), and ABYSS, sort ABYSS \
assemblies, and call blastx and blastn. Only accepts paired end data. Copyright \
2017 by Shawn Rupp, Varsani Lab, Biodesign Institute, Arizona State University. \
This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you \
are welcome to redistribute it under certain conditions.")
	parser.add_argument("-i", help = "Path to manifest file (formatted in \
three tab separated columns: paths to fastqs, batch name \
(all reads for a sample set), sample name (PE reads).")
	parser.add_argument("-o", help = "Path to output directory. All output \
will be written here.")
	parser.add_argument("--noqc", action = "store_true", 
help = "Skip FastQC and Trimmomatic.")
	parser.add_argument("--noblast", action = "store_true",
help = "Ublast/Blast will not be run on sorted contigs.")
	parser.add_argument("--blast", default = False, action = "store_true",
help = "Runs blast on sorted contigs (ublast is run by default).")
	parser.add_argument("--align", action = "store_true",  default = False,
help = "Resume pipeline from blast/ublast.")
	args = parser.parse_args()
	assemble = args.noqc
	if args.align == True:
		assemble = True
		args.noqc = True
	# Load settings and input files
	conf = config(args.noblast)
	fastqs = inputDict(args.i)
	# Make output directory
	outpath = args.o
	if outpath[-1] != "/":
		outpath += "/"
	if not os.path.isdir(outpath):
		os.mkdir(outpath)
	if args.noqc == False:
		# Perform quality control steps
		quals = fastQC(fastqs, outpath, conf["cpu"])
		if quals:
			trim = parseQuality(quals, outpath)
		if trim == True:
			settings = minLen(conf["trim"], conf["k"])
			assemble = trimmomatic(fastqs, outpath, conf["cpu"], settings)
		elif trim == False:
			assemble = True
	if assemble == True:
		if args.align == False:
			contigs = abyss(fastqs, outpath, conf["cpu"], conf["k"], startdir)
		elif args.align == True:
			# Assemble contig paths
			contigs = []
			for batch in fastqs:
				outdir = outpath + "abyss-" + batch + "/"
				contigs.append(outdir + batch + "-8.fa")
	if contigs:
		queries = sort(contigs, outpath, args.align)
	if queries:
		if args.noblast == False:
			if args.blast == True:
				blast(queries, conf)
			else:
				Ublast(queries, conf)
	print(("\tFinished pipeline.\n\tTotal runtime: {}.\n").format(
										datetime.now()-starttime))

if __name__ == "__main__":
	main()
