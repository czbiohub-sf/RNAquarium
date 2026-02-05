#!/usr/bin/env python3
"""
Script to update Diamond search results with missing taxids using local BlastDBCmd
Uses your existing NR_clustered database for perfect compatibility
"""

import gzip
import subprocess
import time
import os
import glob
from pathlib import Path
import tempfile

def lookup_taxids_blastdbcmd(target_ids, blast_db_path, blastdbcmd_path):
    """
    Look up taxids for target IDs using blastdbcmd with your NR database
    
    Args:
        target_ids: list or set of accession IDs
        blast_db_path: path to your nr_cluster_seq database
        blastdbcmd_path: path to blastdbcmd executable
    
    Returns:
        dict mapping target_id -> taxid (or None if not found)
    """
    target_list = list(target_ids)
    taxid_map = {}
    
    print(f"    Looking up {len(target_list)} unique targets using blastdbcmd...")
    
    # Create temporary file with target IDs
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as temp_file:
        temp_filename = temp_file.name
        for target_id in target_list:
            temp_file.write(f"{target_id}\n")
    
    try:
        # Run blastdbcmd in batch mode
        cmd = [
            blastdbcmd_path,
            '-db', blast_db_path,
            '-entry_batch', temp_filename,
            '-outfmt', '%a %T',  # accession and taxid only
        ]
        
        print(f"    Running: {' '.join(cmd)}")
        start_time = time.time()
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode != 0:
            print(f"    Error running blastdbcmd: {result.stderr}")
            # Fall back to individual queries for problematic entries
            return lookup_taxids_individual(target_list, blast_db_path, blastdbcmd_path)
        
        # Parse results
        found_count = 0
        for line in result.stdout.strip().split('\n'):
            if line.strip():
                parts = line.strip().split()
                if len(parts) >= 2:
                    accession = parts[0]
                    taxid = parts[1]
                    
                    # Validate taxid
                    if taxid.isdigit() and int(taxid) > 0:
                        taxid_map[accession] = taxid
                        found_count += 1
                    else:
                        taxid_map[accession] = None
                elif len(parts) == 1:
                    # Only accession returned, no taxid
                    taxid_map[parts[0]] = None
        
        # Mark any targets not found at all
        for target_id in target_list:
            if target_id not in taxid_map:
                taxid_map[target_id] = None
        
        elapsed = time.time() - start_time
        print(f"    Found {found_count}/{len(target_list)} taxids ({found_count/len(target_list)*100:.1f}%) in {elapsed:.2f} seconds")
        
    except subprocess.TimeoutExpired:
        print("    Batch query timed out, falling back to individual queries...")
        return lookup_taxids_individual(target_list, blast_db_path, blastdbcmd_path)
    except Exception as e:
        print(f"    Error with batch query: {e}")
        return lookup_taxids_individual(target_list, blast_db_path, blastdbcmd_path)
    finally:
        # Clean up temp file
        if os.path.exists(temp_filename):
            os.remove(temp_filename)
    
    return taxid_map

def lookup_taxids_individual(target_list, blast_db_path, blastdbcmd_path):
    """
    Fall back to individual blastdbcmd queries
    """
    taxid_map = {}
    found_count = 0
    
    print(f"    Falling back to individual queries for {len(target_list)} targets...")
    
    for i, target_id in enumerate(target_list):
        if (i + 1) % 1000 == 0:
            print(f"      Progress: {i + 1}/{len(target_list)}")
        
        try:
            cmd = [
                blastdbcmd_path,
                '-db', blast_db_path,
                '-entry', target_id,
                '-outfmt', '%a %T'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
            
            if result.returncode == 0 and result.stdout.strip():
                parts = result.stdout.strip().split()
                if len(parts) >= 2 and parts[1].isdigit() and int(parts[1]) > 0:
                    taxid_map[target_id] = parts[1]
                    found_count += 1
                else:
                    taxid_map[target_id] = None
            else:
                taxid_map[target_id] = None
                
        except Exception:
            taxid_map[target_id] = None
    
    print(f"    Individual queries found {found_count}/{len(target_list)} taxids ({found_count/len(target_list)*100:.1f}%)")
    return taxid_map

def process_diamond_file(input_file, blast_db_path, blastdbcmd_path):
    """
    Process a single diamond results file using blastdbcmd lookups
    """
    print(f"Processing {input_file}...")
    
    # Step 1: Read file and collect all rows with N/A taxids
    na_rows = []
    all_lines = []
    unique_targets = set()
    
    print("  Step 1: Reading file and identifying N/A rows...")
    
    # Handle both .gz and .txt files
    if input_file.endswith('.gz'):
        file_opener = lambda f: gzip.open(f, 'rt')
    else:
        file_opener = lambda f: open(f, 'r')
    
    with file_opener(input_file) as f:
        for line_num, line in enumerate(f):
            line = line.strip()
            all_lines.append(line)
            
            if not line:
                continue
                
            fields = line.split('\t')
            if len(fields) < 4:
                continue
                
            # Check if taxname (4th column, index 3) is N/A
            if len(fields) >= 4 and fields[3] == 'N/A':
                target_id = fields[1]  # 2nd column, index 1
                unique_targets.add(target_id)
                na_rows.append((line_num, fields))
    
    if not unique_targets:
        print(f"  No N/A taxids found in {input_file}")
        return
    
    print(f"  Found {len(na_rows)} rows with N/A taxids")
    print(f"  Found {len(unique_targets)} unique target IDs")
    
    # Step 2: Look up taxids using blastdbcmd
    print("  Step 2: Looking up taxids using blastdbcmd...")
    start_time = time.time()
    taxid_map = lookup_taxids_blastdbcmd(unique_targets, blast_db_path, blastdbcmd_path)
    lookup_time = time.time() - start_time
    print(f"    Lookup completed in {lookup_time:.1f} seconds")
    
    # Step 3: Update all N/A rows with the fetched taxids
    print("  Step 3: Updating rows with new taxids...")
    updates_made = 0
    rows_remain_na = 0
    
    for line_idx, fields in na_rows:
        target_id = fields[1]
        new_taxid = taxid_map.get(target_id)
        
        if new_taxid:
            # Update taxid column (3rd column, index 2)
            while len(fields) < 4:
                fields.append('N/A')
            fields[2] = new_taxid
            all_lines[line_idx] = '\t'.join(fields)
            updates_made += 1
        else:
            # No taxid found, ensure proper formatting
            while len(fields) < 4:
                fields.append('N/A')
            if fields[2] in ['', ' ']:
                fields[2] = 'N/A'
            all_lines[line_idx] = '\t'.join(fields)
            rows_remain_na += 1
    
    print(f"  Updated {updates_made}/{len(na_rows)} N/A rows with taxids")
    print(f"  Success rate: {updates_made/len(na_rows)*100:.1f}%")
    
    # Step 4: Write updated file
    print("  Step 4: Writing updated file...")
    
    input_path = Path(input_file)
    if input_path.suffix == '.gz' and input_path.stem.endswith('.txt'):
        base_name = input_path.stem[:-4]  # Remove .txt
        output_file = input_path.parent / f"{base_name}_updated.txt.gz"
    else:
        output_file = input_path.parent / f"{input_path.stem}_updated{input_path.suffix}"
    
    temp_file = str(output_file) + '.tmp'
    
    try:
        with gzip.open(temp_file, 'wt') as f:
            for line in all_lines:
                f.write(line + '\n')
        
        os.replace(temp_file, output_file)
        print(f"  ✓ Updated file saved: {output_file}")
        
    except Exception as e:
        print(f"  Error writing {output_file}: {e}")
        if os.path.exists(temp_file):
            os.remove(temp_file)

def main():
    """
    Main function to process diamond files using blastdbcmd
    """
    # Configuration - UPDATE THESE PATHS
    blast_db_path = "/hpc/scratch/group.data.science/eric_temp/databases/2025/nr_clustered/nr_cluster_seq"  # Full path to your NR database
    blastdbcmd_path = "/hpc/scratch/group.data.science/eric_temp/ncbi-blast-2.15.0+/bin/blastdbcmd"
    
    # Verify paths exist
    if not os.path.exists(f"{blast_db_path}.pal"):
        print(f"Error: BLAST database not found at {blast_db_path}")
        print("Make sure the database path is correct")
        return
    
    if not os.path.exists(blastdbcmd_path):
        print(f"Error: blastdbcmd not found at {blastdbcmd_path}")
        print("Make sure the blastdbcmd path is correct")
        return
    
    # Find files to process
    # FOR TESTING: Single file
    #pattern = "chunk_107_nonzfhum.b.diamond.txt.gz"
    
    # FOR FULL RUN: Uncomment this line
    pattern = "*nonzfhum*.diamond.txt.gz"
    
    files = glob.glob(pattern)
    
    if not files:
        print(f"No files found matching pattern: {pattern}")
        return
    
    print(f"Found {len(files)} files to process")
    print(f"Using BLAST database: {blast_db_path}")
    print(f"Using blastdbcmd: {blastdbcmd_path}")
    
    # Process files
    start_time = time.time()
    for i, file_path in enumerate(sorted(files), 1):
        print(f"\n{'='*60}")
        print(f"Processing file {i}/{len(files)}: {os.path.basename(file_path)}")
        print(f"{'='*60}")
        
        try:
            process_diamond_file(file_path, blast_db_path, blastdbcmd_path)
        except KeyboardInterrupt:
            print("\nInterrupted by user")
            break
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            continue
    
    elapsed = time.time() - start_time
    print(f"\nAll files completed in {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")

if __name__ == "__main__":
    main()
