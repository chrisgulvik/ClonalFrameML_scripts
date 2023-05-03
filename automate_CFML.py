#!/usr/bin/env python
from __future__ import division  #enforce floating point division

__author__ = 'Christopher A. Gulvik'
__version__ = '1.1.0'

import argparse
import fileinput
import getpass
import logging
import os
import re
import socket
import subprocess
import sys
from Bio import AlignIO
from decimal import Decimal, ROUND_UP


def parseArgs():
	parser = argparse.ArgumentParser(
		description='Estimates GTR parameters for ClonalFrameML using PhyML, \
		then infers recombination sites with ClonalFrameML.')
	parser.add_argument('-i', '--fasta', required=True, help='Specify input FASTA file')
	parser.add_argument('-p', '--phyl', required=False, help='Specify input PHYLIP file. IDs must be identical to those in FASTA input. Default generates this file.')
	parser.add_argument('-o', '--outpath', required=False, default='./auto-CFML_output', help='Name of output folder. Default is "./auto-CFML_output".')
	parser.add_argument('-b', '--boots', required=False, type=int, default='100', help='Number of bootstraps for PhyML. Default is "100".')
	parser.add_argument('-s', '--esims', required=False, type=int, default='20', help='Number of uncertainty estimate simulations for CFML. Default is "20".')
	parser.add_argument('-x', '--xmfa', required=False, action='store_true', default=False, help='Input FASTA is XMFA format. Default is False.')
	args = parser.parse_args()
	return args

def dependency(dep):
	""" checks for binary and script availability """
	for path in os.environ.get('PATH', '').split(':'):
		if os.path.exists(os.path.join(path, dep)) and \
		not os.path.isdir(os.path.join(path, dep)):
			return os.path.join(path, dep)
	return None

def fa2phy(fasta, outphy):
	""" converts fasta to phylip """
	with open(fasta) as handle:
		records = AlignIO.parse(handle, 'fasta')
		with open(outphy, 'w') as out:
			AlignIO.write(records, out, 'phylip-relaxed')

def runPhyML(phyl, base, boots, outpath, freqs, kappa, alpha, search, xopts):
	""" runs PhyML to estimate parameters for ClonalFrameML """
	os.system("touch %s/PhyML_%s.log" % (outpath, base))
	os.system("phyml -i %s --datatype nt --bootstrap %s --model GTR -f %s --ts/tv %s --alpha %s --search %s %s | tee %s/PhyML_%s.log" % (phyl, boots, freqs, kappa, alpha, search, xopts, outpath, base))

def getBrDisp(base, ext):
	""" calculates branch dispersion in a newick tree """
	phymltreefile = '%s%s_phyml_tree.txt' % (base, ext)
	with open(phymltreefile) as newick:
		tree = newick.readline()
	l = re.findall(r':[0-9]{1,}\.[E\-0-9]{1,}[\,\(\)]', tree)  # sci notation or simple decimal
	logging.info ('Counted %s leaves in PhyML tree' % len(l))

	branches = []
	for i in l:
		branches.append(float(i.strip(':,)')))
	brDisp = (sum(branches) / len(branches))
	logging.info ('Calculated %s branch dispersion from PhyML tree' % brDisp)
	return brDisp

def getPhyMLest(path, base, ext):
	""" parses PhyML output for variables used as priors for ClonalFrameML """
	phymlstatsfile = '%s%s_phyml_stats.txt' % (base, ext)
	with open(os.path.join(path, phymlstatsfile)) as pstats:
		phymlstats = pstats.read().replace('\n', '')

	# Nucleotides frequencies
	fA = re.search(r'f\(A\)=\s+(0.\d+)', phymlstats)
	fT = re.search(r'f\(T\)=\s+(0.\d+)', phymlstats)
	fC = re.search(r'f\(C\)=\s+(0.\d+)', phymlstats)
	fG = re.search(r'f\(G\)=\s+(0.\d+)', phymlstats)
	freqs = '%s,%s,%s,%s' % (fA.group(1), fC.group(1), fG.group(1), fT.group(1))  # Order (A,C,G,T) matters for PhyML
	logging.info ('Calculated %s frequencies of A,T,C,G from PhyML tree' % freqs)

	# Transition mutations relative to transversion events
	A2C = re.search(r'A <\-> C\s+(\d+.\d+)', phymlstats)
	A2G = re.search(r'A <\-> G\s+(\d+.\d+)', phymlstats)
	A2T = re.search(r'A <\-> T\s+(\d+.\d+)', phymlstats)
	C2G = re.search(r'C <\-> G\s+(\d+.\d+)', phymlstats)
	C2T = re.search(r'C <\-> T\s+(\d+.\d+)', phymlstats)
	G2T = re.search(r'G <\-> T\s+(\d+.\d+)', phymlstats)
	Transitions = (Decimal(A2G.group(1)) + Decimal(C2T.group(1))) # A<->G && C<->T
	logging.info ('Calculated %s transitions from PhyML statistics' % Transitions)
	Transversions = (Decimal(A2C.group(1)) + Decimal(A2T.group(1)) + Decimal(C2G.group(1)) + Decimal(G2T.group(1)))  # All (4) others
	logging.info ('Calculated %s transversions from PhyML statistics' % Transversions)
	TsTv = (Transitions / Transversions).quantize(Decimal('.000001'), rounding=ROUND_UP)

	# Shape of gamma distribution
	gammaShape = re.search(r'Gamma shape parameter:\s+(\d+.\d+)', phymlstats)
	alpha = gammaShape.group(1)
	logging.info ('Calculated alpha as %s from PhyML statistics' % alpha)
	return alpha, TsTv, freqs

def runClonalFrameMLBW(path, base, ext, fasta, outpath, brDisp, kappa, esims):
	""" runs ClonalFrameML with Baum-Welch expectation maximization """
	treefile = '%s%s_phyml_tree.txt' % (base, ext)
	for line in fileinput.input(treefile, inplace=1, backup='.bak'):
		line = re.sub(r'\)[0-9.]*\:', r'):', line.rstrip())
		print line
	cmd = "ClonalFrameML %s%s_phyml_tree.txt %s CFML_%s.BW.out -embranch_dispersion %s -kappa %s -emsim %s | tee %s/CFML_%s.BW.log" % (base, ext, fasta, base, brDisp, kappa, esims, outpath, base)
	subprocess.call(cmd, shell=True)

def getCFMLest(path, base):
	""" parses ClonalFrameML output """
	CFMLstatsfile = 'CFML_%s.BW.out.em.txt' % base
	with open(os.path.join(path, CFMLstatsfile)) as cstats:
		CFMLstats = cstats.read().replace('\n', '')
	print CFMLstats
	# Estimate values on recombination-to-mutation per site, 
	# average length of recombined fragments, and divergence rate of recombination sites
	RperTheta = re.search(r'R\/theta\t(\d+.\d+)', CFMLstats)
	invRecombLen = re.search(r'1\/delta\t(\d+.\d+)', CFMLstats)
	nu = re.search(r'nu\t(\d+.\d+)', CFMLstats)
	logging.info ('Calculated R/Theta as %s from ClonalFrameML statistics' % RperTheta.group(1))
	logging.info ('Calculated 1/RecombLen as %s from ClonalFrameML statistics' % invRecombLen.group(1))
	logging.info ('Calculated nu as %s from ClonalFrameML statistics' % nu.group(1))

	CFMLvals = '"%s %s %s"' % (RperTheta.group(1), invRecombLen.group(1), nu.group(1))
	return CFMLvals

def runClonalFrameMLEM(path, base, ext, fasta, outpath, brDisp, kappa, CFMLvals):
	""" runs ClonalFrameML to estimate recombation parameters for each branch in a tree """
	cmd = "ClonalFrameML %s%s_phyml_tree.txt %s CFML_%s.EM.out -embranch_dispersion %s -kappa %s -embranch true -initial_values %s | tee %s/CFML_%s.EM.log" % (base, ext, fasta, base, brDisp, kappa, CFMLvals, outpath, base)
	subprocess.call(cmd, shell=True)

def main():
	opts = parseArgs()
	
	try: 
		os.makedirs(opts.outpath)
	except OSError:
		if not os.path.isdir(opts.outpath):
			raise

	logging.basicConfig(filename='%s/automate_CFML.log' % opts.outpath, format='%(asctime)s: %(levelname)s: %(message)s', datefmt='%d-%m-%Y %I:%M:%S %p', level=logging.INFO)
	logging.info('%s script initiated by %s on %s with arguments: %s\n' % (sys.argv[0], getpass.getuser(), socket.gethostname(), ' '.join(sys.argv[1:])))

	# Check dependencies
	CFMLbin = dependency('ClonalFrameML')
	PHYMLbin = dependency('phyml')
	if CFMLbin is not None:
		logging.info('found %s' % CFMLbin)
	else:
		print '\tClonalFrameML not found'
		sys.exit(1)
	if PHYMLbin is not None:
		logging.info('found %s' % PHYMLbin)
	else:
		print '\tPhyML not found'
		sys.exit(1)

	path = os.path.dirname(os.path.abspath(opts.fasta))  #handles relational input
	base = os.path.splitext(os.path.basename(opts.fasta))[0]

	# Handle whether a phylip file was also included or not
	if opts.phyl:
		# grabbing EXT enables any phylip extension to work properly (e.g., phy, phyl)
		ext = os.path.splitext(os.path.basename(opts.phyl))[1]
	else:
		fa2phy(opts.fasta, os.path.join(path, base + '.phy'))
		ext = '.phy'
	phyl = os.path.join(path, base + ext)

	if opts.xmfa:
		fasta = '-xmfa_file ' + os.path.realpath(os.path.expanduser(opts.fasta))
	else:
		fasta = os.path.realpath(os.path.expanduser(opts.fasta))

	# The meat
	runPhyML(phyl, base, opts.boots, opts.outpath, "m", "e", "e", "NNI", "-o lr")
	brDisp = getBrDisp(base, ext)
	(alpha, kappa, freqs) = getPhyMLest(path, base, ext)
	runPhyML(phyl, base, opts.boots, opts.outpath, freqs, kappa, alpha, "NNI", " ")
	runClonalFrameMLBW(path, base, ext, fasta, opts.outpath, brDisp, kappa, opts.esims)
	CFMLvals = getCFMLest(path, base)
	runClonalFrameMLEM(path, base, ext, fasta, opts.outpath, brDisp, kappa, CFMLvals)
	logging.info('completed')


if __name__ == '__main__':
	main()

