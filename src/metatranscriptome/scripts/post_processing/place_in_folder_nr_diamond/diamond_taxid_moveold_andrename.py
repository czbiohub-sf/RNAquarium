#!/usr/bin/env python3
"""
Script to organize updated diamond files after taxid update
Moves *_updated.txt.gz files to efetch_update/ subfolder and renames them
"""

import os
import shutil
import glob
from pathlib import Path

def organize_updated_files():
    """
    Find all *_updated.txt.gz files, create efetch_update folder,
    and move files there with original names
    """
    
    # Find all updated files
    updated_files = glob.glob("*_updated.txt.gz")
    
    if not updated_files:
        print("No *_updated.txt.gz files found in current directory")
        return
    
    print(f"Found {len(updated_files)} updated files to organize:")
    for f in updated_files:
        print(f"  {f}")
    
    # Create efetch_update directory if it doesn't exist
    efetch_dir = Path("efetch_update")
    efetch_dir.mkdir(exist_ok=True)
    print(f"\nCreated/verified directory: {efetch_dir}")
    
    # Process each file with safe two-step move+rename
    moved_count = 0
    for updated_file in updated_files:
        try:
            # Generate original filename by removing '_updated' 
            original_name = updated_file.replace('_updated.txt.gz', '.txt.gz')
            
            # Handle the specific pattern like chunk_107_nonzfhum.b.diamond_updated.txt.gz
            if '_updated.txt.gz' in updated_file:
                # Remove _updated part: chunk_107_nonzfhum.b.diamond_updated.txt.gz -> chunk_107_nonzfhum.b.diamond.txt.gz
                original_name = updated_file.replace('_updated.txt.gz', '.txt.gz')
            else:
                print(f"Warning: Unexpected filename pattern: {updated_file}")
                continue
            
            # Step 1: Move file to destination directory with same name
            temp_dest = efetch_dir / updated_file
            shutil.move(updated_file, temp_dest)
            print(f"  Step 1 - Moved: {updated_file} -> {temp_dest}")
            
            # Step 2: Rename to original name within the destination directory
            final_dest = efetch_dir / original_name
            temp_dest.rename(final_dest)
            print(f"  Step 2 - Renamed: {temp_dest.name} -> {final_dest.name}")
            
            moved_count += 1
            
        except Exception as e:
            print(f"Error processing {updated_file}: {e}")
            # If step 1 succeeded but step 2 failed, the file is still safely moved
            # but with the _updated name
    
    print(f"\nSuccessfully organized {moved_count}/{len(updated_files)} files")
    print(f"Files are now in efetch_update/ with original .diamond.txt.gz names")
    
    # Show final directory structure
    if moved_count > 0:
        print(f"\nContents of efetch_update/:")
        for f in sorted(efetch_dir.glob("*.txt.gz")):
            print(f"  {f.name}")

if __name__ == "__main__":
    organize_updated_files()