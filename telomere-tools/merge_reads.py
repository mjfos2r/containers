#!/usr/local/bin/python3

import sys
from Bio import SeqIO

def merge_fastq(filtered_fastq, clipped_fastq, output_fastq):
    """
    Merge a filtered FASTQ file with a clipped FASTQ file.

    Args:
        filtered_fastq (str): FASTQ file of already-filtered reads (unclipped)
        clipped_fastq (str): FASTQ file of clipped reads
        output_fastq (str): Path to write the merged output
    """
    total_filtered = 0
    total_clipped = 0

    with open(output_fastq, 'w') as out_handle:
        # Write clipped reads
        with open(clipped_fastq, 'r') as clipped_handle:
            clipped_records = list(SeqIO.parse(clipped_handle, "fastq"))
            total_clipped = len(clipped_records)
            SeqIO.write(clipped_records, out_handle, "fastq")

        # Write filtered reads
        with open(filtered_fastq, 'r') as filtered_handle:
            filtered_records = list(SeqIO.parse(filtered_handle, "fastq"))
            total_filtered = len(filtered_records)
            SeqIO.write(filtered_records, out_handle, "fastq")

    print(f"Clipped reads: {total_clipped}")
    print(f"Filtered reads: {total_filtered}")
    print(f"Total reads in output: {total_clipped + total_filtered}")

def print_usage():
    print("Usage: python merge_reads.py <filtered_fastq> <clipped_fastq> <output_fastq>")
    print("  filtered_fastq: FASTQ file of pre-filtered reads")
    print("  clipped_fastq: FASTQ file with clipped reads")
    print("  output_fastq: Output merged FASTQ file")
    sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print_usage()

    filtered_fastq = sys.argv[1]
    clipped_fastq = sys.argv[2]
    output_fastq = sys.argv[3]

    merge_fastq(filtered_fastq, clipped_fastq, output_fastq)
