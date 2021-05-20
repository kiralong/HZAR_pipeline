#!/usr/bin/env python3

# Input file
sumstats = 'NW_021939396.1.RDH10.sumstats.tsv'
# Output File
outfile =  './NW_021939396.1.RDH10.parental_fixed_snps.tsv'
out = open(outfile, 'w')
out.write('#Locus_ID\tSNP_column\tChromosome\tBasepair\n')
# Parental population IDs to extract
parent1 = '02SS'
parent2 = '10CG'
# Dictionary to store the parental snp pairs
parental_snps_dict = dict()

# Open and loop over the sumstats file, extract the parental lines and store them as pairs
for line in open(sumstats, 'r'):
    if line[0] == '#':
        continue
    fields = line.strip('\n').split('\t')
    # Extract only the columns we want
    locus_id    = int(fields[0])
    chromosome  = fields[1]
    basepair    = int(fields[2])
    snp_column  = int(fields[3])
    population  = fields[4]
    allele_freq = float(fields[8])
    # Filter lines that are not from the parents
    if population not in [parent1, parent2]:
        continue
    selected_columns = [locus_id, chromosome, basepair, snp_column, population, allele_freq]
    # Create a new SNP ID for the dictionary
    snp_id = '{}_{}'.format(locus_id, snp_column)
    # Load into the dictionary
    parental_snps_dict.setdefault(snp_id, [None, None])
    if population == parent1:
        parental_snps_dict[snp_id][0] = selected_columns
    elif population == parent2:
        parental_snps_dict[snp_id][1] = selected_columns

# Loop over the snps in the dictionary and keep those differentially fixed
for snp in parental_snps_dict:
    snp_pair = parental_snps_dict[snp]
    parent1_snp = snp_pair[0]
    parent2_snp = snp_pair[1]
    # Skip if Parent 1 snp freq is not 0 or 1
    if 0.0 < parent1_snp[5] < 1.0:
        continue
    # Skip if Parent 2  snp freq is not 0 or 1
    if 0.0 < parent2_snp[5] < 1.0:
        continue
    # Compare the pairs for fixes differences
    # Case1: Parent1 = 1, Parent2 = 0
    if parent1_snp[5] == 1.0 and parent2_snp[5] == 0.0:
        row = '{}\t{}\t{}\t{}\n'.format(parent1_snp[0], parent1_snp[3], parent1_snp[1], parent1_snp[2])
        out.write(row)
    # Case2: Parent1 = 0,  Parent2 = 1
    if parent1_snp[5] == 0.0 and parent2_snp[5] == 1.0:
        row = '{}\t{}\t{}\t{}\n'.format(parent1_snp[0], parent1_snp[3], parent1_snp[1], parent1_snp[2])
        out.write(row)

# Close output file
out.close()
