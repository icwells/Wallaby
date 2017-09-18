'''This script will order metagenomic sequencing contigs by length and 
change headers to "sampleID_contig#_lt_seqlength"'''

from collections import OrderedDict
import argparse
import os
import re

def printContigs(outpath, fasta, circ):
	# Prints sorted fasta to output files
	print("\tSaving sorted contigs...")
	with open(outpath + "_min250bp.fasta", "w") as min250:
		with open(outpath + "_all_contigs_sort.fasta", "w") as al:
			for i in fasta:
				# Write remaining sequences
				head = i.replace("~", "_")
				al.write(">" + head + "\n" + fasta[i] + "\n")
				if len(fasta[i]) >= 250:
					min250.write(">" + head + "\n" + fasta[i] + "\n")
		for i in circ:
			if len(circ[i]) >= 250:
				min250.write(i.replace("~", "_") + circ[i] + "\n")
	print("\tFinished.\n")	

def slicer(seq, n):
	# Subset first n nuclotides (potential overlapping end)
	sstring = seq[0:n]
	x = seq.split(sstring)
	L = len(x)
	# Append first n with remaining sequence and last character
	beg = sstring.join(x[:(L-1)])
	end = sstring + x[L-1]
	return beg, end, L

def circularSeqs(fasta, outpath):
	# Identify cicular sequences
	print("\tIdentifying circular contigs...")
	n = 20
	circ = {}
	# Open output fasta and statistic file
	with open(outpath + "_CIRCULAR.fasta", "w") as outfile:
		with open(outpath + "_STAT.txt", "w") as outstat:
			outstat.write("seqid\tnb_hit\tSlen\toverlap\tk\n")
			for i in fasta:
				# Split first n nucleotides and determine lengths of fragments
				beg, end, L3 = slicer(fasta[i], n)
				L1 = len(beg)
				L2 = len(end)
				if L3 > 2 and beg[:L2] == end:
					# Mark circular contigs, write to fasta and stat file
					newid = (">{}_cir_{}\n").format(i, L1).replace("~", "_")
					outfile.write(newid + beg + "\n")
					outstat.write((i +"\t{}\t{}\t{}\t{}\n").format(L3, L1, L2, n))
					circ[newid] = beg
				else:
					# Write stats for all other sequences
					outstat.write((i +"\t{}\t{}\t0\t{}\n").format(L3, L2, n))
	return circ

def sortFasta(infile, filename):
	# Saves input fasta as a sorted dictionary
	print("\n\tSorting contigs...")
	seq = ""
	save = False
	fa = {}
	fasta = OrderedDict()
	# Temporarily remove dashes and underscores
	filename = filename.replace("_", "~")
	filename = filename.replace("-", "~")
	with open(infile, "r") as infile:
		for line in infile:
			line = line.strip()
			if line[0] == ">":
				if save == True:
					fa[head] = seq
					seq = ""
				# Save header without spaces
				head = filename + "_" + line[1:].replace(" ", "_")
				save = True
			else:
				# Append sequence to dict with header as ID
				seq += line
	# Sorts by longest sequence length first; save length in sample ID
	n = 1
	for i in sorted(fa, key=lambda i: len(fa[i]), reverse=True):
		# SampleID_contig#_lt_contigLength
		newid = i.split("_")[0] + "_" + str(n) + "_lt_" +  str(len(fa[i]))
		fasta[newid] = fa[i]
		n += 1		
	return fasta

def main():
	parser = argparse.ArgumentParser(description = 	"filterMetagenomicSequences \
Version 4 (08/28/17) will order contigs by length, update headers, and identify \
circular sequences. Copyright 2017 by Shawn Rupp, Varsani Lab, Biodesign Institute, \
Arizona State University. This program comes with ABSOLUTELY NO WARRANTY. \
This is free software, and you are welcome to redistribute it under certain \
conditions.")
	parser.add_argument("i", help = "Path to input fasta file.")
	args = parser.parse_args()
	infile = args.i
	# Isolate filename
	f = infile[infile.rfind("/")+1:]
	if "-" in f:
		# Remove stage number from ABySS assemblies
		splt = re.split(r"-\d\.", f)
		filename = splt[0]
	else:
		filename = filename[:filename.rfind(".")]
	# Determine output name
	outdir = infile[:infile.rfind("/")+1] + filename + "Sorted/"
	outpath = outdir + filename
	if not os.path.isdir(outdir):
		os.mkdir(outdir)
	fasta = sortFasta(args.i, filename)
	circ = circularSeqs(fasta, outpath)
	printContigs(outpath, fasta, circ)

if __name__ == "__main__":
	main()
