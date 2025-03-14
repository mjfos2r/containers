#!/usr/bin/env python3

import sys
import gzip
from Bio import SeqIO

def extract_original_ids(clipped_fastq):
    """
    Extract the original read IDs from a fastq file with clipped reads.
    These are the IDs of reads that were successfully clipped and should be removed
    from the original dataset.

    Args:
        clipped_fastq (str): Path to the fastq file with clipped reads

    Returns:
        set: Set of original read IDs that were clipped
    """
    original_ids = set()

    # Open the clipped reads file
    with open(clipped_fastq, 'r') as handle:
        for record in SeqIO.parse(handle, "fastq"):
            # Remove the _1 or _2 suffix to get the original read ID
            if record.id.endswith("_1") or record.id.endswith("_2"):
                original_id = record.id[:-2]  # Remove last 2 characters (_1 or _2)
                original_ids.add(original_id)

    return original_ids

def merge_and_filter(original_fastq, clipped_fastq, output_fastq, clipped_ids):
    """
    Merge clipped reads with unclipped reads, removing pre-clipped versions.

    Args:
        original_fastq (str): Path to the original fastq.gz file
        clipped_fastq (str): Path to the clipped reads fastq file
        output_fastq (str): Path to write the merged output
        clipped_ids (set): Set of original read IDs that were clipped
    """
    # Counters for stats
    total_reads = 0
    removed_reads = 0
    kept_unclipped = 0

    # First write all clipped reads to the output file
    clipped_records = []
    with open(clipped_fastq, 'r') as handle:
        clipped_records = list(SeqIO.parse(handle, "fastq"))

    with open(output_fastq, 'w') as out_handle:
        # Write all clipped reads to the output file
        SeqIO.write(clipped_records, out_handle, "fastq")

        # Now process the original fastq and keep only the unclipped reads
        with gzip.open(original_fastq, 'rt') as in_handle:
            for record in SeqIO.parse(in_handle, "fastq"):
                total_reads += 1

                # Skip reads that were clipped (we've already added their clipped versions)
                if record.id in clipped_ids:
                    removed_reads += 1
                else:
                    # Keep reads that weren't clipped
                    SeqIO.write(record, out_handle, "fastq")
                    kept_unclipped += 1

    print(f"Total original reads processed: {total_reads}")
    print(f"Unclipped reads kept: {kept_unclipped}")
    print(f"Clipped reads: {len(clipped_records)}")
    print(f"Reads Removed (unclipped): {removed_reads}")
    print(f"Total reads in output: {kept_unclipped + len(clipped_records)}")

def print_usage():
    print("Usage: python merge_clipped_reads.py <original_fastq.gz> <clipped_fastq> <output_fastq>")
    print("  original_fastq.gz: Original FASTQ.gz file with all reads")
    print("  clipped_fastq: FASTQ file with the clipped reads (from previous script)")
    print("  output_fastq: Output FASTQ file containing both clipped fragments and unclipped reads")
    sys.exit(1)

if __name__ == "__main__":
    # Check if correct number of arguments is provided
    if len(sys.argv) != 4:
        print_usage()

    # Get command line arguments
    original_fastq = sys.argv[1]  # Original gzipped fastq file
    clipped_fastq = sys.argv[2]   # Fastq file with clipped reads from previous script
    output_fastq = sys.argv[3]    # Output file for merged reads

    # Extract the IDs of reads that were successfully clipped
    clipped_ids = extract_original_ids(clipped_fastq)
    print(f"Found {len(clipped_ids)} unique reads that were successfully clipped")

    # Merge clipped reads with unclipped reads, removing pre-clipped versions
    merge_and_filter(original_fastq, clipped_fastq, output_fastq, clipped_ids)