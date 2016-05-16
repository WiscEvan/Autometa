#!/usr/bin/env python

# Program to assess the quality of deconvolution of a simulated or defined metagenome dataset
# Inputs:
# 	1. Sam file of reads aligned to the known reference genomes
#	2. Table of species classifications of each contig in the reference genome fasta
#	3. Table of read ranges for each species
#	4. Sam file of reads aligned to the deconvoluted assembly
#	5. Table of bin classifications of each contig in the assembly

# Outputs:
#	1. 'Binning accuracy' table, header: bin\tgenome\tpercent
#		Tells you the percent of each bin that belongs to different genomes
#	2. 'Binning recovery' table, header: genome\tbin\tpercent
#		Tells you the percent of each genome found in a bin
#	3. 'Chimera' table, header: contig\tgenome\tpercent
#		Tells you the percent of each contig that belongs to different genomes

import sys
import getopt
import gzip

def is_alignment_congruent_with_ref(read_name, contig_aligned_to, read_ranges, contig_species):
	species = contig_species[contig_aligned_to]
	if does_read_belong(read_name, species, read_ranges):
		return True
	else:
		return False

def does_read_belong(read_name, species, read_ranges):
	# Read name assumed to be in the form read_235
	read_list = read_name.split('_')
	read_number = read_list[1]
	species_start_read = read_ranges[species]['start']
	species_end_read = read_ranges[species]['end']
	if read_number >= species_start_read and read_number <= species_end_read:
			return True
	else:
		return False

def species_classification_for_read(read_name, read_ranges, excluded_reads):
	classified_species = None
	for species in read_ranges:
		if does_read_belong(read_name, species, read_ranges):
			if read_name not in excluded_reads:
				classified_species = species
	return classified_species

def get_species_percents(read_counts):
	percents = {}
	total_read_count = 0
	for species in read_counts:
		total_read_count += read_counts[species]
	if total_read_count == 0:
		print 'get_species_percents: total_read_count = 0'
		sys.exit(2)
	for species in read_counts:
		percent = (read_counts[species]/total_read_count)*100
		percents[species] = percent
	return percents

ref_sam_path = None # -r --refsam
ref_species_table_path = None # -s --refspecies, table made by make_simulated_metagenome_control_fasta.py
ref_read_ranges_table_path = None # -q --readranges, table made by make_simulated_metagenome.py
asm_sam_path = None # -a --asmsam
bin_classifications_table_path = None # -b --bintable
output_prefix = None
bin_column = None

opts,args = getopt.getopt(sys.argv[1:],"hr:s:q:a:b:o:c",["help", "refsam=", "refspecies=", "readranges=", "asmsam=", "bintable=", "outputprefix=", "column="])

for opt, arg in opts:
	if opt in ('-h', '--help'):
		print 'assess_metagenome_deconvolution.py -r <sam alignment to ref genomes> -s <ref contig species table> -q <ref read table> -a <sam alignment to assembly> -b <bin classification table for assembly>'
		sys.exit()
	elif opt in ('-r', '--refsam'):
		ref_sam_path = arg
	elif opt in ('-s', '--refspecies'):
		ref_species_table_path = arg
	elif opt in ('-q', '--readranges'):
		ref_read_ranges_table_path = arg
	elif opt in ('-a', '--asmsam'):
		asm_sam_path = arg
	elif opt in ('-b', '--bintable'):
		bin_classifications_table_path = arg
	elif opt in ('-o', '--outputprefix'):
		output_prefix = arg
	elif opt in ('-c', '--column'):
		bin_column = arg

if not bin_column:
	bin_column = 'db.cluster' # Default value for bin column header to look out for

print 'Reference SAM: ' + ref_sam_path
print 'Contig species table: ' + ref_species_table_path
print 'Read ranges table: ' + ref_read_ranges_table_path
print 'Assembly SAM: ' + asm_sam_path
print 'Bin classifications table: ' + bin_classifications_table_path
print 'Bin column: ' + bin_column
print 'Output prefix: ' + output_prefix
print '\n'

# 1. Parse read ranges table, so that we can spot non-unique reads in the reference alignment
print 'Parsing read ranges table...'
ranges = {} # Dictionary of dictionaries, keyed by species
range_table_rows = ((row.rstrip('\n')) for row in open(ref_read_ranges_table_path))
for i,row in enumerate(range_table_rows):
	if not i == 0:
		row_list = row.split('\t')
		ranges[row_list[0]] = { 'start':row_list[2], 'end':row_list[3] }

# 2. Go through reference contig table, and remember which species each contig belongs to
print 'Parsing contig species table...'
species = {} # Dictionary, keyed by contig, stores species
species_table_rows = ((row.rstrip('\n')) for row in open(ref_species_table_path))
for i,row in enumerate(species_table_rows):
	if not i == 0:
		row_list = row.split('\t')
		species[row_list[0]] = row_list[1]

# 3. Go through reference sam file, and flag non-unique reads
# Non-unique reads are spotted in the sam file as reads that occur in genomes other than their originating genome
# (because they must be in their originating genome, therefore seeing one outside means they occur in at least two)
# Note: this does not flag up reads that occur more than once in their originating genomes
print 'Finding non-unique reads in reference SAM...'
non_unique_reads = {} # Dictionary that just contains reads found in more than one genome

# If sam file is a gz file, use gzip, otherwise normal open
if ref_sam_path[-3:] == '.gz':
	with gzip.open(ref_sam_path, 'rb') as ref_sam:
		for line in ref_sam:
			# Skip header section, where lines begin with '@'
			if not line[0] == '@':
				line_list = line.split('\t')
				print line_list
				read_name = line_list[0]
				contig_name = line_list[1]
				contig_species = species[contig_name]
				# Find if read "belongs" to the species that this contig is from
				if not does_read_belong(read_name, contig_species, ranges):
					non_unique_reads[read_name] = 1
else:
	with open(ref_sam_path) as ref_sam:
		for line in ref_sam:
			# Skip header section, where lines begin with '@'
			if not line[0] == '@':
				line_list = line.split('\t')
				read_name = line_list[0]
				contig_name = line_list[1]
				contig_species = species[contig_name]
				# Find if read "belongs" to the species that this contig is from
				if not does_read_belong(read_name, contig_species, ranges):
					non_unique_reads[read_name] = 1

# 3a. Make a data structure that records the number of unique reads for each bin.
print 'Working out how many unique reads there are per genome...'
number_of_unique_reads = {} # Dictionary keyed by bin name
for genome in ranges:
	start_read = ranges[genome]['start']
	end_read = ranges[genome]['end']
	start_read_list = start_read.split('_')
	end_read_list = start_read.split('_')
	start_read_number = int(start_read_list[1])
	end_read_number = int(end_read_list[1])
	counter = start_read_number
	unique_reads = 0
	while counter <= end_read_number:
		read_name = 'read_' + counter
		if read_name not in non_unique_reads:
			unique_reads += 1
		counter += 1
	number_of_unique_reads[genome] = unique_reads

# 4. Go through assembly sam file, and count read classifications for each contig
print 'Parsing assembly SAM, counting species reads...'
contig_classifications = {} # Dictionary of dictionaries, which will hold running tallies of reads assigned to different species

# If sam file is a gz file, use gzip
if asm_sam_path[-3:] == '.gz':
	with gzip.open(asm_sam_path, 'rb') as asm_sam:
		for line in asm_sam:
			# Skip header section, where lines begin with '@'
			if not line[0] == '@':
				line_list = line.split('\t')
				read_name = line_list[0]
				contig_name = line_list[1]
				read_species = species_classification_for_read(read_name, ranges, non_unique_reads) # Returns None if a non-unique read
				if read_species:
					if contig_name in contig_classifications:
						if read_species in contig_classifications[contig_name]:
							contig_classifications[contig_name][read_species] += 1
						else:
							contig_classifications[contig_name][read_species] = 1
					else:
						contig_classifications[contig_name] = { read_species: 1 }
else:
	with open(asm_sam_path) as asm_sam:
		for line in asm_sam:
			# Skip header section, where lines begin with '@'
			if not line[0] == '@':
				line_list = line.split('\t')
				read_name = line_list[0]
				contig_name = line_list[1]
				read_species = species_classification_for_read(read_name, ranges, non_unique_reads) # Returns None if a non-unique read
				if read_species:
					if contig_name in contig_classifications:
						if read_species in contig_classifications[contig_name]:
							contig_classifications[contig_name][read_species] += 1
						else:
							contig_classifications[contig_name][read_species] = 1
					else:
						contig_classifications[contig_name] = { read_species: 1 }

# 5. We now have enough information to write a table showing how chimeric contigs are
# Output table in the format contig\tgenome\treads\tpercent
chimera_table_path = output_prefix + '_chimera_table'
print 'Writing chimera table ' + chimera_table_path + '...'
chimera_table = open(chimera_table_path, 'w')
chimera_table.write('contig\tgenome\treads\tpercent\n')
for contig in contig_classifications:
	percents = get_species_percents(contig_classifications[contig])
	for species in contig_classifications[contig]:
		chimera_table.write(contig + '\t' + species + '\t' + contig_classifications[contig][species] + '\t' + percents[species] + '\n')
chimera_table.close

# 6. We need to go through the bin table to make a datastructure containing the classification of each contig
print 'Making bin datastructure...'
contig_bins = {} # Dictionary, keyed by contig, stores bin classifications
bin_table_rows = ((row.rstrip('\n')) for row in open(bin_classifications_table_path))

# Find out which column to use for bin classification 
bin_column_index = None
number_found = 0
first_line_list = bin_table_rows[0].split('\t')
for i,value in enumerate(first_line_list):
	if value == bin_column:
		bin_column_index = i
		number_found += 1
if number_found > 1:
	print 'Error, bin table has more than one column headed ' + bin_column
	sys.exit(2)
if not bin_column_index:
	print 'Error, could not find column ' + bin_column + ' in bin table'
	sys.exit(2)

contig_column_index = None
number_found = 0
for i,value in enumerate(first_line_list):
	if value == 'contig':
		contig_column_index = i
		number_found += 1
if number_found > 1:
	print 'Error, bin table has more than one contig column'
	sys.exit(2)
if not contig_column_index:
	print 'Error, could not find contig column in bin table'
	sys.exit(2)

for i,row in enumerate(range_table_rows):
	if not i == 0:
		contig = row[contig_column_index]
		bin_name = row[bin_column_index]
		contig_bins[contig] = bin_name

# Now make a data structure that totals up the species reads for each bin
bin_classifications = {} # Dictionary keyed by bins, then species
for contig in contig_bins:
	current_bin = contig_bins[contig]
	if current_bin not in bin_classifications:
		bin_classifications[current_bin] = {}

	for species in contig_classifications[contig]:
		number_reads = contig_classifications[contig][species]
		if species in bin_classifications[current_bin]:
			bin_classifications[current_bin][species] += number_reads
		else:
			bin_classifications[current_bin][species] = number_reads

# 7. Make 'Binning accuracy' table, header: bin\tgenome\treads\tpercent
bin_accuracy_table_path = output_prefix + '_bin_accuracy_table'
print 'Writing binning accuracy table ' + bin_accuracy_table_path + '...'
bin_accuracy_table = open(bin_accuracy_table_path, 'w')
bin_accuracy_table.write('bin\tgenome\treads\tpercent\n')
for bin_name in bin_classifications:
	percents = get_species_percents(bin_classifications[bin_name])
	for species in bin_classifications[bin_name]:
		bin_accuracy_table.write(bin_name + '\t' + species + '\t' + bin_classifications[bin_name][species] + '\t' + percents[species] + '\n')
bin_accuracy_table.close

# 8. Make a 'Binning recovery' table, header: genome\tbin\treads\tpercent
print 'Counting genome reads in bins...'
genome_reads_in_bins = {}
for bin_name in bin_classifications:
	for species in bin_classifications[bin_name]:
		if species not in genome_reads_in_bins:
			genome_reads_in_bins[species] = {}

		if bin_name in genome_reads_in_bins[species]:
			genome_reads_in_bins[species][bin_name] += bin_classifications[bin_name][species]
		else:
			genome_reads_in_bins[species][bin_name] = bin_classifications[bin_name][species]

bin_recovery_table_path = output_prefix + '_bin_recovery_table'
print 'Writing binning recovery table ' + bin_recovery_table_path + '...'
bin_recovery_table = open(bin_recovery_table_path, 'w')
bin_recovery_table.write('genome\tbin\treads\tpercent\n')
for species in genome_reads_in_bins:
	# Percents are calculated based on the previously calculated number of unique reads per reference genome
	# because - not all the genome might end up in bins/assembled
	percents = {}
	for bin_name in genome_reads_in_bins[species]:
		number_of_reads = genome_reads_in_bins[species][bin_name]
		percent = (number_of_reads / number_of_unique_reads[species])*100
		percents[bin_name] = percent

	for bin_name in genome_reads_in_bins[species]:
		bin_recovery_table.write(species + '\t' + bin_name + '\t' + genome_reads_in_bins[species][bin_name] + '\t' + percents[bin_name] + '\n')
bin_recovery_table.close