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
                    utrS, utrE, grpN, siteType, miRinSpecies, grpTy, \
                    spInGrp, spInGrpWtype, ORF_overlap = l.split('\t')
    
            msaKey = ':'.join([msaS, msaE])            #  Multiple-sequence alignment start and end joined
            posKey = ':'.join([utrS, utrE])            #  Species-specific UTR start and stop coordinates
            tgtKey = ':'.join([geneTrans, miRgroup])   #  Join gene name, transcript name, and miR

            if tgtKey not in tgtSpList:                #  If tgtKey is not in tgtSpList, add it
                tgtSpList[tgtKey] = {}
                
   # If sites perfectly overlap, are contained, or encompass another site,
   # count all sites as the same location
            for msa_key in tgtSpList[tgtKey]:
                kmsaS, kmsaE = msa_key.split(':')      #  Split dictionary msa key into start and stop  
                kmsaS, kmsaE, msaS, msaE = map(int, [kmsaS, kmsaE, msaS, msaE])   # Convert to integer
                if msaS >= kmsaS and kmsaE >= msaE or kmsaS >= msaS and msaE >= kmsaE:
                    msaKey = msa_key
                    tgtSpList[tgtKey][msaKey].append(sID)   # Add species to species list for that gene/location
                    break
            else:
                tgtSpList[tgtKey][msaKey] = [sID]           # Gene/location not in tgtSpList, so add it
    
    # Only record target if in species of interest
            if sID == species:                              # If the species ID == species of interest
                try:
                    tgtUTR2MSA[tgtKey][posKey] = msaKey     # Add the gene/miR/UTR position == gene/msa position dict
                    tgtSite[tgtKey][posKey] = siteType      # Add the gene/miR/UTR position miR site type to dict
                except KeyError:
                    tgtUTR2MSA[tgtKey] = {posKey : msaKey}
                    tgtSite[tgtKey] = {posKey : siteType}
    return tgtUTR2MSA, tgtSite, tgtSpList


def write_output(output, tgtUTR2MSA, tgtSite, tgtSpList):
    '''Write output to file, this is the near completed scorecard (needs UTR len)'''
    print 'Output saved as:', output
    with open(output, 'w') as f:
        for tgtKey in tgtSite:                        # for gene/miR
            gene, trans, miRgroup = tgtKey.split(':') # split into parts
            for pos in tgtSite[tgtKey]:               # for position in gene/mir/UTR position dict
                msaKey = tgtUTR2MSA[tgtKey][pos]      # get the corresponding msa position
                try:
                    spList = ':'.join(tgtSpList[tgtKey][msaKey]) # get the species that share that site
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
