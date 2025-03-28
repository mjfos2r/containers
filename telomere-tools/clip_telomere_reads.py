#!/usr/local/bin/python3
import sys
import gzip
from Bio import SeqIO
import statistics

def get_double_motif_hit_coordinates_dict(bedfile) -> dict:
    reads_dict = {}
    positive_strand_reads = set() # create set for positive strand reads
    negative_strand_reads = set() # create set for negative strand reads
    with open(bedfile, 'r') as file: # open the bedfile
        for line in file: # loop through each line
            # split the line into fields and grab the read_id and strand
            fields = line.strip().split('\t')
            read_id = fields[0]
            strand = fields[5]
            # grab the start and end coordinates of the motif hit
            start = int(fields[1])
            end = int(fields[2])
            length = end - start
            # if the read id is not in the dictionary, create an empty list for that id
            if reads_dict.get(read_id) is None:
                reads_dict[read_id] = []
            # extract the coordinates of the motif hits for each strand
            if strand == '+':
                positive_strand_reads.add(read_id)
                reads_dict[read_id].append({"+" : (start, end)})
            elif strand == '-':
                negative_strand_reads.add(read_id)
                reads_dict[read_id].append({"-" : (start, end)})
    # get the reads that have motif hits on both strands
    paired_reads = positive_strand_reads.intersection(negative_strand_reads)
    # filter the reads dictionary to only include reads with motif hits on both strands
    reads_dict = {k: v for k, v in reads_dict.items() if k in paired_reads}
    return reads_dict

def parseFastqDict(reads_file, reads_dict) -> dict:
    reads = {}
    with gzip.open(reads_file, "rt") as reads_in:
        for record in SeqIO.parse(reads_in, "fastq"): # parse fastq with biopython
            # filter the reads to only include those with motif hits on both strands
            if record.id in reads_dict.keys():
                reads[record.id] = record
            else:
                continue
    return reads

def clipReads(motif_reads, motif_hits) -> list: # list of SeqObjects
    clipped_reads = []
    # loop through the reads with motif hits on both strands
    for read in motif_reads:
        current_read = motif_hits[read]
        # loop through the motif hits for each read
        for hit in current_read:
            # get the coordinates of the motif hits for each strand
            if hit.get("+"):
                positive_hit = hit["+"]
                continue
            elif hit.get("-"):
                neg_hit = hit["-"]
            if abs(max(positive_hit) - min(neg_hit)) > 50: # this is double the 21bp motif length + 8bp
                continue
            else:
                telomeres = positive_hit + neg_hit
                # get cut site based on shared coordinate.
                cut_site = round(statistics.median(telomeres))
                current_record = motif_reads[read]
                # clip the reads at the cut site
                clipped_1 = current_record[cut_site::]
                clipped_2 = current_record[:cut_site:]
                #rename clipped fragments
                clipped_1.id = current_record.id + "_1"
                clipped_2.id = current_record.id + "_2"

                clipped_reads.append(clipped_1)
                clipped_reads.append(clipped_2)
    return clipped_reads

def writeFastq(clipped_reads, output):
    # take list of SeqObjects and write to fastq file
    with open(output, "w") as output_handle:
        SeqIO.write(clipped_reads, output_handle, "fastq")

if __name__ == "__main__":
    # Check if correct number of arguments is provided
    if len(sys.argv) != 4:
        print_usage()

    # Get command line arguments
    bedfile = sys.argv[1]
    reads_file = sys.argv[2]
    output_file = sys.argv[3]

    # get the coordinates of the motif hits for all reads that have adjacent hits on separate strands
    motif_hits = get_double_motif_hit_coordinates_dict(bedfile)
    # parse fastq file and get reads with motif hits
    motif_reads = parseFastqDict(reads, motif_hits)
    # clip reads at the point of symmetry between the two motif hits
    clipped_reads = clipReads(motif_reads, motif_hits)
    # write the clipped reads to a new fastq file
    writeFastq(clipped_reads, output)