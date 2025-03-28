#!/usr/local/bin/python3

import gzip
import sys
from Bio import SeqIO
import statistics

def get_double_motif_hit_coordinates_dict(bedfile) -> dict:
    reads_dict = {}
    with open(bedfile, 'r') as file:
        for line in file:
            fields = line.strip().split('\t')
            read_id = fields[0]
            strand = fields[5]
            start = int(fields[1])
            end = int(fields[2])
            if reads_dict.get(read_id) is None:
                reads_dict[read_id] = []
            reads_dict[read_id].append({strand: (start, end)})

    # Keep reads with at least two hits on the same strand
    filtered_reads = {}
    for read_id, hits in reads_dict.items():
        strand_counts = {'+': 0, '-': 0}
        for hit in hits:
            for strand in hit:
                strand_counts[strand] += 1
        if strand_counts['+'] >= 2 or strand_counts['-'] >= 2:
            filtered_reads[read_id] = hits

    return filtered_reads

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

def clipReads(motif_reads, motif_hits) -> list:
    clipped_reads = []

    for read_id, hits in motif_hits.items():
        current_record = motif_reads[read_id]

        plus_hits = [hit["+"] for hit in hits if "+" in hit]
        minus_hits = [hit["-"] for hit in hits if "-" in hit]

        if len(plus_hits) >= 2:
            same_strand_hits = plus_hits
        elif len(minus_hits) >= 2:
            same_strand_hits = minus_hits
        else:
            continue

        # Sort hits by start coord
        same_strand_hits = sorted(same_strand_hits, key=lambda x: x[0])

        # Use first two hits for now
        hit1 = same_strand_hits[0]
        hit2 = same_strand_hits[1]

        # Determine the gap between the motifs using one end and one start
        trim_coord_1 = hit1[1]  # end of hit1
        trim_coord_2 = hit2[0]  # start of hit2


        # Make sure coord1 < coord2 for correct midpoint calculation
        lower, upper = sorted([trim_coord_1, trim_coord_2])
        if upper - lower > 50: # this realistically shouldn't be over 30-40bp but let's make it 50 for safety
            continue

        # Midpoint between end of one hit and start of the other
        cut_site = round((lower + upper) / 2)

        # Clip the read
        clipped_1 = current_record[cut_site:]
        clipped_2 = current_record[:cut_site]

        # append to read ids
        clipped_1.id = current_record.id + "_1"
        clipped_2.id = current_record.id + "_2"
        clipped_1.description = ""
        clipped_2.description = ""

        clipped_reads.append(clipped_1)
        clipped_reads.append(clipped_2)

    return clipped_reads


def writeFastq(clipped_reads, output):
    # take list of SeqObjects and write to fastq file
    with open(output, "w") as output_handle:
        SeqIO.write(clipped_reads, output_handle, "fastq")

def print_usage():
    print("Usage: python script.py <bedfile> <reads_file> <output_file>")
    print("  bedfile: BED file with motif hits")
    print("  reads_file: FASTQ.gz file with reads")
    print("  output_file: Output FASTQ file for clipped reads")
    sys.exit(1)

if __name__ == "__main__":
    # Check if correct number of arguments is provided
    if len(sys.argv) != 4:
        print_usage()

    # Get command line arguments
    bedfile = sys.argv[1]
    reads_file = sys.argv[2]
    output_file = sys.argv[3]

    # Get the coordinates of the motif hits for all reads that have adjacent hits on separate strands
    motif_hits = get_double_motif_hit_coordinates_dict(bedfile)

    # Parse fastq file and get reads with motif hits
    motif_reads = parseFastqDict(reads_file, motif_hits)

    # Clip reads at the point of symmetry between the two motif hits
    clipped_reads = clipReads(motif_reads, motif_hits)

    # Write the clipped reads to a new fastq file
    writeFastq(clipped_reads, output_file)