#!/usr/bin/env python3
"""
Clean and resolve Autocycler GFA graphs.

Usage:
    python clean_gfa.py <input.gfa> <output_prefix> [options]

Options:
    --min_depth FLOAT       Remove nodes at or below this depth
    --max_depth FLOAT       Remove nodes above this depth
    --split_telos NODES     Comma-separated telomere nodes to split at center
    --remove NODES          Comma-separated nodes to remove
    --keep left|right       Which half to keep when splitting [default: left]
    --min_length INT        Minimum contig length for FASTA output [default: 0]

Example:
    python clean_gfa.py C34_cleaned.gfa C34_final --min_depth 2 --max_depth 24 --split_telos 29,30,74,75 --remove 50,51 --min_length 2500
"""

import sys
import argparse
import subprocess
import shutil
from datetime import datetime


class Logger:
    """Simple logger that writes to both stdout and a log file."""
    def __init__(self, log_path):
        self.log_path = log_path
        self.log_file = open(log_path, 'w')
        self.log(f"=== GFA Cleaning Log ===")
        self.log(f"Started: {datetime.now().isoformat()}")
        self.log("")

    def log(self, message):
        print(message)
        self.log_file.write(message + '\n')
        self.log_file.flush()

    def close(self):
        self.log("")
        self.log(f"Finished: {datetime.now().isoformat()}")
        self.log_file.close()


def parse_depth(fields):
    """Extract depth from DP:f:X or dp:f:X tag."""
    for field in fields:
        if field.upper().startswith('DP:F:'):
            return float(field.split(':')[2])
    return None


def load_gfa(path):
    """Load GFA into segments, links, and other lines."""
    segments = {}
    links = []
    other = []

    with open(path) as f:
        for line in f:
            if line.startswith('S\t'):
                parts = line.strip().split('\t')
                segments[parts[1]] = {'seq': parts[2], 'tags': parts[3:]}
            elif line.startswith('L\t'):
                links.append(line.strip().split('\t'))
            else:
                other.append(line)

    return segments, links, other


def write_gfa(path, segments, links, other):
    """Write GFA from segments, links, and other lines."""
    with open(path, 'w') as out:
        for line in other:
            out.write(line)

        for node_id, data in segments.items():
            tags = '\t'.join(data['tags']) if data['tags'] else ''
            if tags:
                out.write(f"S\t{node_id}\t{data['seq']}\t{tags}\n")
            else:
                out.write(f"S\t{node_id}\t{data['seq']}\n")

        for link in links:
            out.write('\t'.join(link) + '\n')


def dedupe_links(links, logger=None):
    """Remove duplicate links."""
    seen = set()
    unique = []
    dupes = []
    for link in links:
        key = tuple(link[1:5])
        if key not in seen:
            seen.add(key)
            unique.append(link)
        else:
            dupes.append(key)

    if len(links) != len(unique) and logger:
        logger.log(f"Deduplicated links: {len(links)} -> {len(unique)}")
        for d in dupes:
            logger.log(f"  Removed duplicate: {d}")

    return unique


def filter_depth(segments, links, min_depth=None, max_depth=None, logger=None):
    """Remove nodes outside depth range."""
    keep_nodes = set()
    removed_min = []
    removed_max = []

    for node_id, data in segments.items():
        depth = parse_depth(data['tags'])
        seq_len = len(data['seq'])

        keep = True
        if depth is not None:
            if min_depth is not None and depth < min_depth:
                removed_min.append((node_id, seq_len, depth))
                keep = False
            elif max_depth is not None and depth > max_depth:
                removed_max.append((node_id, seq_len, depth))
                keep = False

        if keep:
            keep_nodes.add(node_id)

    if logger:
        if min_depth is not None and removed_min:
            logger.log(f"Min depth filter (depth <= {min_depth}): removed {len(removed_min)} nodes")
            for node_id, seq_len, depth in sorted(removed_min, key=lambda x: -x[2]):
                logger.log(f"  Node {node_id}: {seq_len} bp, {depth}x")

        if max_depth is not None and removed_max:
            logger.log(f"Max depth filter (depth > {max_depth}): removed {len(removed_max)} nodes")
            for node_id, seq_len, depth in sorted(removed_max, key=lambda x: -x[2]):
                logger.log(f"  Node {node_id}: {seq_len} bp, {depth}x")

        logger.log(f"Nodes remaining: {len(keep_nodes)}")

    # Filter segments
    new_segments = {k: v for k, v in segments.items() if k in keep_nodes}

    # Filter links
    new_links = [l for l in links if l[1] in keep_nodes and l[3] in keep_nodes]

    return new_segments, new_links


def split_telomeres(segments, links, telo_nodes, keep='left', logger=None):
    """Split telomere nodes at center, keeping one half."""
    if not telo_nodes:
        return segments, links

    # Validate nodes exist
    missing = set(telo_nodes) - set(segments.keys())
    if missing:
        if logger:
            logger.log(f"Warning: telomere nodes not found: {missing}")
        telo_nodes = [n for n in telo_nodes if n not in missing]

    if not telo_nodes:
        return segments, links

    if logger:
        logger.log(f"Splitting telomere nodes (keeping {keep} half):")

    # Find max node ID
    max_id = max(int(n) for n in segments.keys())

    # Map old node IDs to new
    node_mapping = {}

    for target_node in telo_nodes:
        new_id = str(max_id + 1)
        max_id += 1

        seq = segments[target_node]['seq']
        mid = len(seq) // 2
        kept_seq = seq[:mid] if keep == 'left' else seq[mid:]

        if logger:
            logger.log(f"  Node {target_node} ({len(seq)} bp) -> Node {new_id} ({len(kept_seq)} bp)")
            logger.log(f"    Original: {seq}")
            logger.log(f"    Kept:     {kept_seq}")

        node_mapping[target_node] = {
            'new_id': new_id,
            'seq': kept_seq,
            'tags': segments[target_node]['tags']
        }

    # Build new segments
    new_segments = {}
    for node_id, data in segments.items():
        if node_id in node_mapping:
            info = node_mapping[node_id]
            new_segments[info['new_id']] = {'seq': info['seq'], 'tags': info['tags']}
        else:
            new_segments[node_id] = data

    # Update links
    new_links = []
    for link in links:
        from_node, from_orient, to_node, to_orient = link[1:5]
        overlap = link[5] if len(link) > 5 else '0M'

        skip = False
        new_from = from_node
        new_to = to_node

        if from_node in node_mapping:
            info = node_mapping[from_node]
            if keep == 'left' and from_orient == '+':
                skip = True
            elif keep == 'right' and from_orient == '-':
                skip = True
            else:
                new_from = info['new_id']

        if to_node in node_mapping:
            info = node_mapping[to_node]
            if keep == 'left' and to_orient == '-':
                skip = True
            elif keep == 'right' and to_orient == '+':
                skip = True
            else:
                new_to = info['new_id']

        if not skip:
            new_links.append(['L', new_from, from_orient, new_to, to_orient, overlap])

    return new_segments, new_links


def remove_nodes(segments, links, remove_list, logger=None):
    """Remove specified nodes."""
    if not remove_list:
        return segments, links

    remove_set = set(remove_list)
    missing = remove_set - set(segments.keys())
    if missing and logger:
        logger.log(f"Warning: nodes to remove not found: {missing}")

    removed = remove_set & set(segments.keys())

    if logger and removed:
        logger.log(f"Manually removing {len(removed)} nodes:")
        for node_id in sorted(removed, key=lambda x: int(x) if x.isdigit() else x):
            data = segments[node_id]
            depth = parse_depth(data['tags'])
            depth_str = f", {depth}x" if depth else ""
            logger.log(f"  Node {node_id}: {len(data['seq'])} bp{depth_str}")

    # Filter segments
    new_segments = {k: v for k, v in segments.items() if k not in remove_set}

    # Filter links
    new_links = [l for l in links if l[1] not in remove_set and l[3] not in remove_set]

    return new_segments, new_links


def run_autocycler_clean(in_gfa, out_gfa, duplicate=None, logger=None):
    """Run autocycler clean to merge linear paths."""
    cmd = ['autocycler', 'clean', '-i', in_gfa, '-o', out_gfa]
    if duplicate:
        cmd.extend(['-d', duplicate])
    if logger:
        logger.log(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        if logger:
            logger.log(f"Warning: autocycler clean failed: {result.stderr}")
        return False
    if logger and result.stderr:
        for line in result.stderr.strip().split('\n'):
            logger.log(f"  {line}")
    return True


def run_gfa2fasta(in_gfa, out_fasta, logger=None):
    """Run autocycler gfa2fasta."""
    cmd = ['autocycler', 'gfa2fasta', '-i', in_gfa, '-o', out_fasta]
    if logger:
        logger.log(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        if logger:
            logger.log(f"Warning: autocycler gfa2fasta failed: {result.stderr}")
        return False
    return True


def filter_fasta_by_length(in_fasta, out_fasta, min_length, logger=None):
    """Filter FASTA by minimum length using seqkit."""
    cmd = ['seqkit', 'seq', '-m', str(min_length), in_fasta]
    if logger:
        logger.log(f"Filtering FASTA by min length {min_length} bp")
    with open(out_fasta, 'w') as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        if logger:
            logger.log(f"Warning: seqkit filter failed: {result.stderr}")
        return False
    return True


def get_fasta_stats(fasta_path):
    """Get stats from seqkit."""
    result = subprocess.run(['seqkit', 'stats', '-T', fasta_path],
                          capture_output=True, text=True)
    if result.returncode == 0:
        return result.stdout
    return None


def main():
    parser = argparse.ArgumentParser(description='Clean and resolve Autocycler GFA graphs.')
    parser.add_argument('input_gfa', help='Input GFA file')
    parser.add_argument('output_prefix', help='Output prefix for GFA and FASTA files')
    parser.add_argument('--min_depth', type=float, help='Remove nodes at or below this depth')
    parser.add_argument('--max_depth', type=float, help='Remove nodes above this depth')
    parser.add_argument('--split_telos', help='Comma-separated telomere nodes to split')
    parser.add_argument('--remove', help='Comma-separated nodes to remove')
    parser.add_argument('--keep', choices=['left', 'right'], default='left',
                        help='Which half to keep when splitting telomeres [default: left]')
    parser.add_argument('--duplicate', help='Comma-separated nodes to duplicate for TIR resolution.')
    parser.add_argument('--min_length', type=int, default=0,
                        help='Minimum contig length for FASTA output [default: 0]')
    args = parser.parse_args()

    # Initialize logger
    log_path = f"{args.output_prefix}.log"
    logger = Logger(log_path)

    # Log parameters
    logger.log("Parameters:")
    logger.log(f"  Input: {args.input_gfa}")
    logger.log(f"  Output prefix: {args.output_prefix}")
    if args.min_depth is not None:
        logger.log(f"  Min depth: {args.min_depth}")
    if args.max_depth is not None:
        logger.log(f"  Max depth: {args.max_depth}")
    if args.split_telos:
        logger.log(f"  Split telomeres: {args.split_telos}")
    if args.remove:
        logger.log(f"  Remove nodes: {args.remove}")
    logger.log(f"  Keep half: {args.keep}")
    if args.min_length > 0:
        logger.log(f"  Min contig length: {args.min_length}")
    logger.log("")

   # Parse node lists
    telo_nodes = args.split_telos.replace(' ', '').split(',') if args.split_telos else []
    remove_list = args.remove.replace(' ', '').split(',') if args.remove else []

    logger.log(f"Loading {args.input_gfa}")
    segments, links, other = load_gfa(args.input_gfa)
    logger.log(f"Loaded {len(segments)} segments, {len(links)} links")
    logger.log("")

    # Dedupe links first
    links = dedupe_links(links, logger)

    # Step 1: Filter by MAX depth (removes tangles, preserves node IDs for manual ops)
    if args.max_depth is not None:
        logger.log("")
        segments, links = filter_depth(segments, links, min_depth=None, max_depth=args.max_depth, logger=logger)

    # Step 2: Split telomeres (user-specified node IDs)
    if telo_nodes:
        logger.log("")
        segments, links = split_telomeres(segments, links, telo_nodes, args.keep, logger)

    # Step 3: Remove specified nodes (user-specified node IDs)
    if remove_list:
        logger.log("")
        segments, links = remove_nodes(segments, links, remove_list, logger)

    # Dedupe links again after modifications
    logger.log("")
    links = dedupe_links(links, logger)

    # Write intermediate GFA
    intermediate_gfa = f"{args.output_prefix}_intermediate.gfa"
    write_gfa(intermediate_gfa, segments, links, other)
    logger.log(f"Wrote intermediate GFA: {intermediate_gfa}")

    # Step 4: Run autocycler clean to merge paths
    logger.log("")
    cleaned_gfa = f"{args.output_prefix}_cleaned.gfa"
    if run_autocycler_clean(intermediate_gfa, cleaned_gfa, duplicate=args.duplicate, logger=logger):
        logger.log(f"Wrote cleaned GFA: {cleaned_gfa}")
    else:
        logger.log("Skipping autocycler clean, using intermediate")
        shutil.copy(intermediate_gfa, cleaned_gfa)

    # Step 5: Filter by MIN depth (after autocycler clean, node IDs may have changed)
    final_gfa = f"{args.output_prefix}.gfa"
    if args.min_depth is not None:
        logger.log("")
        segments, links, other = load_gfa(cleaned_gfa)
        segments, links = filter_depth(segments, links, min_depth=args.min_depth, max_depth=None, logger=logger)
        links = dedupe_links(links, logger)
        write_gfa(final_gfa, segments, links, other)
        logger.log(f"Wrote final GFA: {final_gfa}")
    else:
        shutil.copy(cleaned_gfa, final_gfa)
        logger.log(f"Wrote final GFA: {final_gfa}")

    # Step 6: Convert to FASTA
    logger.log("")
    temp_fasta = f"{args.output_prefix}_temp.fasta"
    if run_gfa2fasta(final_gfa, temp_fasta, logger):
        # Step 7: Filter by length if specified
        final_fasta = f"{args.output_prefix}.fasta"
        if args.min_length > 0:
            filter_fasta_by_length(temp_fasta, final_fasta, args.min_length, logger)
            logger.log(f"Wrote final FASTA: {final_fasta}")
        else:
            shutil.move(temp_fasta, final_fasta)
            logger.log(f"Wrote final FASTA: {final_fasta}")

    # Stats
    logger.log("")
    logger.log("=== Final Assembly Stats ===")
    stats = get_fasta_stats(f"{args.output_prefix}.fasta")
    if stats:
        for line in stats.strip().split('\n'):
            logger.log(line)

    logger.log("")
    logger.log(f"Log written to: {log_path}")
    logger.close()


if __name__ == '__main__':
    main()
