#!/usr/bin/python2

###############################################################################
#
#  Usage:
#    python addUTRlength.py input UTRfile [-sp species -fo output]
#
#  Where:
#    input = parsed targetscan file
#    UTRfile = species UTR file
#    species = species ID number (Default = 9606)
#    output = output file name (Default = species_scorecard.txt)
#
###############################################################################

import argparse

parser = argparse.ArgumentParser(description='''Parse TargetScan output.''', 
        formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('INPUT', 
        action='store', 
        help='Input file. This is the output from parse_target_scan.py')
parser.add_argument('UTR_FILE', 
        action='store', 
        help='UTR file which was generated earlier in pipeline')
parser.add_argument('-fo', 
        action='store', 
        dest='fo', 
        default='species_scorecard.txt',
        help='Name for the output file')
parser.add_argument('-sp', 
        action='store', 
        dest='species', 
        default=9606, 
        help='Species ID for the species of interest')
arg = parser.parse_args()

tgtFile = arg.INPUT
utrFile = arg.UTR_FILE
outFile = arg.fo
species = arg.species

print 'Input file is:', tgtFile
print 'UTR input file is:', utrFile
print 'Species of interest is:', species
print 'Output file is:', outFile

LenH, genLen, genT = {}, {}, {}
with open(utrFile, 'r') as f:
    for l in f:
        genTrans, spID, seq = l.rstrip().split('\t')
        gene, trans = genTrans.split(':')
        if spID == species:
            seqNoAlign = seq.replace('-', '')
            length = len(seqNoAlign)
            LenH[trans] = length
            if gene not in genLen or length > genLen[gene]:
                genLen[gene] = length
                genT[gene] = trans

with open(tgtFile, 'r') as f, open(outFile, 'w') as fo:
    for l in f:
        gene, trans, miRgroup, siteType, pos, consList = l.rstrip().split('\t')
        if gene not in genT:
            print "Error: {} should be in genT, {}".format(gene, genT[gene])
        if genT[gene] == trans:
            try:
                utrLen = str(LenH[trans])
            except KeyError:
                print "Error: missing", trans
            line = '\t'.join([gene.upper(), trans, utrLen, miRgroup, siteType, pos, consList])
            fo.write('{}\n'.format(line)) 
