#!/usr/bin/python2

###############################################################################
#
#  Parse the target scan output into the correct format for the scorecard.
#  This includes isolating only those mRNAs in which the species of interest is
#  bound by the miRNA.
#
#  Usage:
#     python parse_target_scan.py INPUT [-fo OUTPUT -sp SPECIES]
#  Where:
#     INPUT = target scan output
#     OUTPUT = output file name (Default = parsed_target_scan.txt)
#     SPECIES = ID of species of interest (Default = 9606)
#
###############################################################################

import argparse


def parse_targetscan_by_species_ID(input, species):
    '''Parse the target scan output lines and filter by species ID'''
    print 'Input file is:', input
    print 'Species of interest is:', species

    tgtSite, tgtSpList, tgtUTR2MSA = {}, {}, {}
    
    with open(input, 'r') as f:
        header = f.readline().rstrip()
        for l in f:
    # Parse target scan output line
            geneTrans, miRgroup, sID, msaS, msaE, \
                    utrS, utrE, grpN, siteTy, miRinSpecies, grpTy, \
                    spInGrp, spInGrpWtype, ORF_overlap = l.split('\t')
    
            msaKey = ':'.join([msaS, msaE])
            posKey = ':'.join([utrS, utrE])
            tgtKey = ':'.join([geneTrans, miRgroup])

            if tgtKey not in tgtSpList:
                tgtSpList[tgtKey] = {}
                
   # If sites perfectly overlap, are contained, or encompass another site,
   # count all sites as the same location
            for msa_key in tgtSpList[tgtKey]:
                kmsaS, kmsaE = msa_key.split(':')
                kmsaS, kmsaE, msaS, msaE = map(int, [kmsaS, kmsaE, msaS, msaE])
                if msaS >= kmsaS and kmsaE >= msaE or kmsaS >= msaS and msaE >= kmsaE:
                    msaKey = msa_key
                    tgtSpList[tgtKey][msaKey].append(sID)
                    break
            else:
                tgtSpList[tgtKey][msaKey] = [sID]
    
    # Only record target if in species of interest
            if sID == species:
                try:
                    tgtUTR2MSA[tgtKey][posKey] = msaKey
                    tgtSite[tgtKey][posKey] = siteTy
                except KeyError:
                    tgtUTR2MSA[tgtKey] = {posKey : msaKey}
                    tgtSite[tgtKey] = {posKey : siteTy}
    return tgtUTR2MSA, tgtSite, tgtSpList


def write_output(output, tgtUTR2MSA, tgtSite, tgtSpList):
    '''Write output to file, this is the near completed scorecard (needs UTR len)'''
    print 'Output saved as:', output
    with open(output, 'w') as f:
        for tgtKey in tgtSite:
            gene, trans, miRgroup = tgtKey.split(':')
            for pos in tgtSite[tgtKey]:
                msaKey = tgtUTR2MSA[tgtKey][pos]
                try:
                    spList = ':'.join(tgtSpList[tgtKey][msaKey])
                except KeyError:
                    print tgtKey, pos
    
                output_line = '\t'.join([gene, trans, miRgroup, tgtSite[tgtKey][pos], pos, spList])
                f.write('{}\n'.format(output_line))


def main(arg):
    UTR2MSA, Site, SpList = parse_targetscan_by_species_ID(arg.INPUT, arg.species)
    write_output(arg.fo, UTR2MSA, Site, SpList)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Parse TargetScan output.''', 
            formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('INPUT', 
            action='store', 
            help='Input file. This is the output from TargetScan')
    parser.add_argument('-fo', 
            action='store', 
            dest='fo', 
            default='parsed_target_scan.txt',
            help='Name for the output file')
    parser.add_argument('-sp', 
            action='store', 
            dest='species', 
            default=9606, 
            help='Species ID for the species of interest')
    main(parser.parse_args())
