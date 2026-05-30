#!/bin/bash

#Script to add taxid_lca_NTorNR column to large TSV file
#Handles various special character cases with exact and inexact searches & Uses single AWK approach - loads taxonomy once, processes everything in memory - Fastest approach for large datasets

# Configuration
TSV_FILE="/path/to/RNAquarium_outputs/taxonomy_hits_viruses0_fullcols_mostrecent.tsv.gz"
NAMES_DMP="/path/to/databases/taxonomizr/aug2025taxonomy/names.dmp"
OUTPUT_FILE="/path/to/RNAquarium_outputs/taxonomy_hits_viruses0_fullcols_mostrecent_withtaxids.tsv.gz"
TEMP_DIR="temp_taxid_processing_viruses0"

echo "Processing large TSV file: $TSV_FILE"
echo "Using taxonomy database: $NAMES_DMP"
echo "Output file: $OUTPUT_FILE"

# Validate input files
if [ ! -f "$TSV_FILE" ]; then
    echo "Error: TSV file not found: $TSV_FILE"
    exit 1
fi

if [ ! -f "$NAMES_DMP" ]; then
    echo "Error: Taxonomy file not found: $NAMES_DMP"
    exit 1
fi

# Create temporary directory
mkdir -p "$TEMP_DIR"
cd "$TEMP_DIR"

# Step 1: Extract unique taxnames from column 7
echo "Step 1: Preparing single AWK processing..."
echo "This approach loads all taxonomy data into memory once for maximum speed"


echo "Step 1: Extracting unique taxnames from column 7..."
echo "Processing 30M+ rows - this may take 5-10 minutes..."

zcat "$TSV_FILE" | \
  tail -n +2 | \
  cut -f7 | \
  grep -v '^$' | \
  sort -u > unique_taxnames.txt

UNIQUE_COUNT=$(wc -l < unique_taxnames.txt)
echo "Found $UNIQUE_COUNT unique taxnames"



# Step 3
# Categorize for statistics (not used in processing, just for reporting)
echo "Analyzing taxname categories for statistics..."


# 1. Underscore cases (highest priority)
grep '_' unique_taxnames.txt > underscore_cases.txt

# 2. Bracket cases - check remaining names
grep -v '_' unique_taxnames.txt | grep '\[' > bracket_cases.txt

# 3. Asterisk cases - check names without underscores or brackets  
grep -v '_' unique_taxnames.txt | grep -v '\[' | grep '\*' > asterisk_cases.txt

# 4. General inexact cases - remaining problematic characters
grep -v '_' unique_taxnames.txt | grep -v '\[' | grep -v '\*' | \
  grep '[()\/\\^$?|-]' > inexact_cases.txt

# 5. Exact search cases - no problematic characters
#grep -v '[_*\[\]()\/\\^$?|-]' unique_taxnames.txt > exact_search.txt
# 5. Exact search cases - exclude ALL problematic characters (fix the exclusion chain)
grep -v '_' unique_taxnames.txt | grep -v '\[' | grep -v '\*' | \
  grep -v '[()\/\\^$?|-]' > exact_search.txt

echo "Category analysis:"
echo "  Exact search candidates: $(wc -l < exact_search.txt)"
echo "  Underscore cases: $(wc -l < underscore_cases.txt)"
echo "  Bracket cases: $(wc -l < bracket_cases.txt)"
echo "  Asterisk cases: $(wc -l < asterisk_cases.txt)"
echo "  General inexact cases: $(wc -l < inexact_cases.txt)"

echo ""
echo "Step 3: Creating comprehensive AWK script..."
echo "Loading entire taxonomy database into memory for fast lookups..."

# Create comprehensive AWK script for taxonomy processing
cat > process_taxonomy.awk << 'EOF'
BEGIN {
    FS = OFS = "\t"
    print "Loading taxonomy database into memory..." > "/dev/stderr"
}

# First AWK section: Load taxonomy database (names.dmp)
FNR == NR {
    # Parse names.dmp format: taxid |tab| |tab| name |tab| |tab| name_class
    split($0, fields, "\t\\|\t")
    if (length(fields) >= 2) {
        taxname = fields[2]
        taxid = fields[1]
        
        # Store exact matches (case-sensitive) - log multiple matches
        if (!exact_tax[taxname]) {
            exact_tax[taxname] = taxid
        } else if (exact_tax[taxname] != taxid) {
            print "Multiple exact matches for: " taxname " (using " exact_tax[taxname] ", also found " taxid ")" > "/dev/stderr"
        }
        
        # Store inexact matches (case-insensitive) - log multiple matches
        lowercase_name = tolower(taxname)
        if (!inexact_tax[lowercase_name]) {
            inexact_tax[lowercase_name] = taxid
        } else if (inexact_tax[lowercase_name] != taxid) {
            print "Multiple inexact matches for: " lowercase_name " (using " inexact_tax[lowercase_name] ", also found " taxid ")" > "/dev/stderr"
        }
    }
    next
}

# Second AWK section: Process TSV data
FNR < NR {
    # Header row
    if (FNR == 1) {
        # Add new column header after column 7
        output_line = ""
        for (i = 1; i <= 7; i++) {
            output_line = output_line $i (i < 7 ? OFS : "")
        }
        output_line = output_line OFS "taxid_lca_NTorNR"
        for (i = 8; i <= NF; i++) {
            output_line = output_line OFS $i
        }
        print output_line
        next
    }
    
    # Data rows - validate minimum columns
    if (NF < 7) {
        print "Warning: Row " FNR " has only " NF " columns, skipping" > "/dev/stderr"
        next
    }
    
    # Get taxname from column 7
    original_taxname = $7
    found_taxid = ""
    
    # FIRST PASS: Try different lookup strategies with precedence order
    
    # 1. Exact match first
    if (original_taxname in exact_tax) {
        found_taxid = exact_tax[original_taxname]
    }
    # 2. Handle underscore cases (highest priority for special chars)
    else if (index(original_taxname, "_") > 0) {
        modified_name = original_taxname
        # Remove words containing underscore and preceding spaces
        gsub(/ [^ ]*_[^ ]*/, "", modified_name)
        gsub(/^[^ ]*_[^ ]*/, "", modified_name)
        # Clean up spaces
        gsub(/^ +| +$/, "", modified_name)
        
        if (modified_name != "" && length(modified_name) > 0) {
            lowercase_modified = tolower(modified_name)
            if (lowercase_modified in inexact_tax) {
                found_taxid = inexact_tax[lowercase_modified]
            }
        }
    }
    # 3. Handle bracket cases - use escaping for inexact search
    else if (index(original_taxname, "[") > 0) {
        # Escape brackets and do inexact search on original name
        escaped_name = original_taxname
        gsub(/\[/, "\\[", escaped_name)
        gsub(/\]/, "\\]", escaped_name)
        lowercase_escaped = tolower(escaped_name)
        
        if (lowercase_escaped in inexact_tax) {
            found_taxid = inexact_tax[lowercase_escaped]
        }
    }
    # 4. Handle asterisk cases - use escaping for inexact search
    else if (index(original_taxname, "*") > 0) {
        # Escape asterisks and do inexact search on original name
        escaped_name = original_taxname
        gsub(/\*/, "\\*", escaped_name)
        lowercase_escaped = tolower(escaped_name)
        
        if (lowercase_escaped in inexact_tax) {
            found_taxid = inexact_tax[lowercase_escaped]
        }
    }
    # No step 5, no final else
#     # 5. General inexact search for problematic characters (UPDATED: includes dashes)
#     else if (match(original_taxname, /[()\/\\^$?|-]/)) {
#         lowercase_original = tolower(original_taxname)
#         if (inexact_tax[lowercase_original]) {
#             found_taxid = inexact_tax[lowercase_original]
#         }
#     }
    # 6. Clean names - inexact search (REMOVE the else, make it conditional)
    # 6. Remove this section entirely - let everything else go to second pass

    
    # SECOND PASS: If still not found, try broader matching strategies
    # SECOND PASS: If still not found, try simple inexact search first
    if (found_taxid == "") {
        print "DEBUG: Second pass for: " original_taxname > "/dev/stderr"
        # ... rest of second pass
        lowercase_original = tolower(original_taxname)
        split(lowercase_original, words, " ")
        # Most important: try unchanged taxname with inexact lookup
        if (lowercase_original in inexact_tax) {
            print "DEBUG: Found simple inexact match for: " original_taxname > "/dev/stderr"
            found_taxid = inexact_tax[lowercase_original]
        }
        
        # Step 2: Try first 3 words as exact search (fast)
        if (found_taxid == "" && length(words) >= 3) {
            three_word_name = words[1] " " words[2] " " words[3]
            if (three_word_name in inexact_tax) {
                print "DEBUG: Found 3-word match for: " original_taxname > "/dev/stderr"
                found_taxid = inexact_tax[three_word_name]
            }
        }
        
        # Step 3: Try first 2 words as exact search (fast)
        if (found_taxid == "" && length(words) >= 2) {
            two_word_name = words[1] " " words[2]
            if (two_word_name in inexact_tax) {
                print "DEBUG: Found 2-word match for: " original_taxname > "/dev/stderr"
                found_taxid = inexact_tax[two_word_name]
            }
        }
        
#         # Final fallback: genus-only prefix matching
#         if (found_taxid == "" && length(words) >= 2) {
#             # Try "genus species" prefix matching
#             two_word_name = words[1] " " words[2]
#             # Look for taxonomy entries that START with this pattern
#             for (tax_name in inexact_tax) {
#                 if (index(tax_name, two_word_name) == 1) {  # Starts with genus_species
#                     found_taxid = inexact_tax[tax_name]
#                     print "DEBUG: Found prefix match for: " original_taxname " -> " found_taxid " (matched: " tax_name ")" > "/dev/stderr"
#                     break
#                 }
#             }
#         }
    }
    
    # Use found taxid or N/A
    final_taxid = (found_taxid != "") ? found_taxid : "N/A"
    
    # Build output line with new column after column 7
    output_line = ""
    for (i = 1; i <= 7; i++) {
        output_line = output_line $i (i < 7 ? OFS : "")
    }
    output_line = output_line OFS final_taxid
    for (i = 8; i <= NF; i++) {
        output_line = output_line OFS $i
    }
    print output_line
}
EOF


# After processing, add logging section
echo ""
echo "Step 4: Checking for multiple match cases..."
echo "Multiple taxid matches are logged during taxonomy database loading"
echo "Check stderr output above for any 'Multiple matches' warnings"

# Execute the AWK script with logging
echo "Executing comprehensive taxonomy processing with second pass fallback..."
echo "This may take 10-30 minutes depending on file size..."

# Capture both stdout and stderr to see multiple match warnings
# Process: names.dmp first (loaded into memory), then TSV data (processed)
awk -f process_taxonomy.awk "$NAMES_DMP" <(zcat "$TSV_FILE") 2> processing_warnings.log | gzip > "$OUTPUT_FILE"

echo "AWK processing completed!"

# Step 4: Generate detailed statistics
echo "Step 4b: Generating statistics..."

# Extract just the new taxid column to analyze success rate
echo "Analyzing success rate from processed data..."
zcat "$OUTPUT_FILE" | tail -n +2 | cut -f8 > new_taxid_column.txt

TOTAL_ROWS=$(wc -l < new_taxid_column.txt)
FOUND_COUNT=$(grep -v '^N/A$' new_taxid_column.txt | wc -l)
NOT_FOUND_COUNT=$(grep '^N/A$' new_taxid_column.txt | wc -l)

if [ "$TOTAL_ROWS" -gt 0 ]; then
    SUCCESS_RATE=$(echo "scale=1; $FOUND_COUNT * 100 / $TOTAL_ROWS" | bc -l 2>/dev/null || echo "unknown")
else
    SUCCESS_RATE="0"
fi

echo "Part I processing statistics:"
echo "  Total data rows processed: $TOTAL_ROWS"
echo "  Taxids found: $FOUND_COUNT"
echo "  Taxids not found (N/A): $NOT_FOUND_COUNT"
echo "  Success rate: $SUCCESS_RATE%"

# Show breakdown by category (based on original analysis)
echo ""
echo "Expected category breakdown (from preprocessing):"
echo "  Unique taxnames analyzed: $UNIQUE_COUNT"
echo "  Exact search candidates: $(wc -l < exact_search.txt 2>/dev/null || echo 0)"
echo "  Underscore cases: $(wc -l < underscore_cases.txt 2>/dev/null || echo 0)"
echo "  Bracket cases: $(wc -l < bracket_cases.txt 2>/dev/null || echo 0)"
echo "  Asterisk cases: $(wc -l < asterisk_cases.txt 2>/dev/null || echo 0)"
echo "  General inexact cases: $(wc -l < inexact_cases.txt 2>/dev/null || echo 0)"

# Show examples of matches and non-matches
echo ""
echo "Sample successful assignments:"
zcat "$OUTPUT_FILE" | tail -n +2 | awk -F'\t' '$8 != "N/A" {print "  " $7 " -> " $8}' | head -5

echo ""
echo "Sample failed assignments:"
zcat "$OUTPUT_FILE" | tail -n +2 | awk -F'\t' '$8 == "N/A" {print "  " $7 " -> N/A"}' | head -5


# Show multiple match statistics
if [ -s processing_warnings.log ]; then
# Count unique taxnames with multiple matches
    UNIQUE_MULTIPLE_EXACT=$(grep "Multiple exact matches" processing_warnings.log | cut -d':' -f2 | cut -d'(' -f1 | sort -u | wc -l)
    UNIQUE_MULTIPLE_INEXACT=$(grep "Multiple inexact matches" processing_warnings.log | cut -d':' -f2 | cut -d'(' -f1 | sort -u | wc -l)
    MULTIPLE_EXACT=$(grep "Multiple exact matches" processing_warnings.log | wc -l)
    MULTIPLE_INEXACT=$(grep "Multiple inexact matches" processing_warnings.log | wc -l)
    cd ..
    cp "$TEMP_DIR/processing_warnings.log" "multiple_matches_and_partI_assignments.log"
    echo ""
    echo "Multiple match summary:"
    echo "  Multiple exact matches found: $UNIQUE_MULTIPLE_EXACT"
    echo "  Multiple inexact matches found: $UNIQUE_MULTIPLE_INEXACT"
    echo "  Details saved to: multiple_matches_and_partI_assignments.log"
else
    echo "No multiple matches detected during processing"
fi

# Extract everything up to the opening parenthesis
grep "Multiple exact matches for:" multiple_matches_and_partI_assignments.log | \
  sed 's/Multiple exact matches for: //' | \
  sed 's/ (using.*$//' | \
  sort -u > log_exact_taxnames.txt

grep "Multiple inexact matches for:" multiple_matches_and_partI_assignments.log | \
  sed 's/Multiple inexact matches for: //' | \
  sed 's/ (using.*$//' | \
  sort -u > log_inexact_taxnames.txt

# Combine and find intersection with your TSV taxnames
cat log_exact_taxnames.txt log_inexact_taxnames.txt | \
  sort -u > log_all_taxnames.txt

# Find which ones are actually in your TSV data
comm -12 <(sort unique_taxnames.txt) log_all_taxnames.txt > tsv_multiple_candidates.txt

echo "TSV taxnames with multiple possible matches: $(wc -l < tsv_multiple_candidates.txt)"

## rename outputs
mv log_exact_taxnames.txt log_exact_taxnames_viruses.txt
mv log_inexact_taxnames.txt log_inexact_taxnames_viruses.txt
mv log_all_taxnames.txt log_all_taxnames_viruses.txt
mv tsv_multiple_candidates.txt tsv_multiple_candidates_viruses.txt

# Step 5: Cleanup and verification
#cd ..
echo ""
echo "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

echo ""
echo "Part I Processing completed successfully!"
echo "Original file: $TSV_FILE"  
echo "New file: $OUTPUT_FILE"

# Show file sizes
ORIGINAL_SIZE=$(du -h "$TSV_FILE" | cut -f1 2>/dev/null || echo "unknown")
NEW_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1 2>/dev/null || echo "unknown")
echo "File size: $ORIGINAL_SIZE -> $NEW_SIZE"

# Verify new column structure
echo ""
echo "Verifying column structure (first 3 rows):"
echo "Columns 6-9 (should show: col6, taxname_lca_NTorNR, taxid_lca_NTorNR, col9):"
zcat "$OUTPUT_FILE" | head -3 | cut -f6-9

# Final validation
echo ""
echo "Final Part I validation:"
HEADER_COLS=$(zcat "$OUTPUT_FILE" | head -1 | awk -F'\t' '{print NF}')
DATA_COLS=$(zcat "$OUTPUT_FILE" | head -2 | tail -1 | awk -F'\t' '{print NF}')
echo "  Header columns: $HEADER_COLS"
echo "  Data columns: $DATA_COLS"

if [ "$HEADER_COLS" -eq "$DATA_COLS" ]; then
    echo "  ✓ Column count validation passed"
else
    echo "  ⚠ Warning: Column count mismatch"
fi

# Check if new column exists
NEW_COL_HEADER=$(zcat "$OUTPUT_FILE" | head -1 | cut -f8)
if [ "$NEW_COL_HEADER" = "taxid_lca_NTorNR" ]; then
    echo "  ✓ New column header found in position 8"
else
    echo "  ⚠ Warning: Expected 'taxid_lca_NTorNR' in column 8, found: '$NEW_COL_HEADER'"
fi

echo ""
echo "Complete! The taxid_lca_NTorNR column has been added after taxname_lca_NTorNR."
echo "Performance: Single AWK pass - loads taxonomy once, processes entire dataset in memory."
echo ""
echo "PART I of Script completed successfully!"

echo ""
echo "============================================================"
echo "PART II: Processing remaining N/A cases using tax_species_NTorNR"
echo "============================================================"

## note for nonhost.tsv.gz tax_species_NTorNR is in column 53 but in viruses0_fullcols it is in column 107!! ACTUALLY 108 AFTER ADDING TAXID COLUMN!!
# Define Part II output file with proper extension handling for full paths
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
BASE_NAME=$(basename "$OUTPUT_FILE" .tsv.gz)
PART2_OUTPUT_FILE="${OUTPUT_DIR}/${BASE_NAME}_final.tsv.gz"

# Create new temp directory for Part II processing
TEMP_DIR_PART2="temp_taxid_processing_part2"
mkdir -p "$TEMP_DIR_PART2"
cd "$TEMP_DIR_PART2"

echo "Part II input file: $OUTPUT_FILE"
echo "Part II output file: $PART2_OUTPUT_FILE"

# Step 1: Validate column structure and extract N/A rows
echo "Part II Step 1: Validating file structure and extracting N/A rows..."

# Check if output file exists and has content
if [ ! -f "$OUTPUT_FILE" ]; then
    echo "Error: Part I output file not found: $OUTPUT_FILE"
    cd ..
    rm -rf "$TEMP_DIR_PART2"
    exit 1
fi

# Validate column count
HEADER_COLS=$(zcat "$OUTPUT_FILE" | head -1 | awk -F'\t' '{print NF}')
if [ "$HEADER_COLS" -lt 108 ]; then
    echo "Error: File has only $HEADER_COLS columns, need at least 108 for tax_species_NTorNR"
    cd ..
    rm -rf "$TEMP_DIR_PART2"
    exit 1
fi

echo "File validation passed: $HEADER_COLS columns found"

# Extract rows with N/A in column 8
zcat "$OUTPUT_FILE" | awk -F'\t' 'NR==1 {print; next} $8=="N/A" {print}' > na_rows.tsv


# Count N/A rows (excluding header)
NA_ROWS=$(tail -n +2 na_rows.tsv | wc -l)
echo "Found $NA_ROWS rows with N/A taxids"

if [ "$NA_ROWS" -eq 0 ]; then
    echo "No N/A rows found - Part II not needed"
    # Copy original file as final output
    cp "$OUTPUT_FILE" "$PART2_OUTPUT_FILE"
    cd ..
    rm -rf "$TEMP_DIR_PART2"
    echo "Final output: $PART2_OUTPUT_FILE (identical to Part I output)"
else
    # Step 2: Extract unique species names from column 108
    echo "Part II Step 2: Extracting unique species names from column 108..."
    tail -n +2 na_rows.tsv | cut -f108 | grep -v '^$' | grep -v '^N/A$' | sort -u > unique_species_names.txt
    
    UNIQUE_SPECIES_COUNT=$(wc -l < unique_species_names.txt)
    echo "Found $UNIQUE_SPECIES_COUNT unique species names for lookup"
    
    if [ "$UNIQUE_SPECIES_COUNT" -eq 0 ]; then
        echo "No valid species names found in column 108 - copying original file"
        cp "$OUTPUT_FILE" "$PART2_OUTPUT_FILE"
        cd ..
        rm -rf "$TEMP_DIR_PART2"
    else
        # Step 3: Create species name to taxid mapping
        echo "Part II Step 3: Creating species name lookup mapping..."
        
        # Validate names.dmp accessibility
        if [ ! -f "$NAMES_DMP" ]; then
            echo "Error: Cannot access taxonomy database: $NAMES_DMP"
            cd ..
            rm -rf "$TEMP_DIR_PART2"
            exit 1
        fi
        
        # Create AWK script for species-based lookups
        cat > species_lookup.awk << 'EOF'
        
BEGIN {
    FS = OFS = "\t"
    print "Loading taxonomy database for species lookup..." > "/dev/stderr"
}

# Load taxonomy database
FNR == NR {
    split($0, fields, "\t\\|\t")
    if (length(fields) >= 2) {
        taxname = fields[2]
        taxid = fields[1]
        
        # Store lowercase for inexact matching - take first occurrence only
        lowercase_name = tolower(taxname)
        if (!species_tax[lowercase_name]) {
            species_tax[lowercase_name] = taxid
        }
    }
    next
}

# Process species names file
FNR < NR {
    species_name = $1
    lowercase_species = tolower(species_name)
    found_species_taxid = ""
    
    # Try exact inexact lookup
    if (lowercase_species in species_tax) {
        found_species_taxid = species_tax[lowercase_species]
    }
    
    # Output mapping: original_name -> taxid
    final_species_taxid = (found_species_taxid != "") ? found_species_taxid : "N/A"
    print species_name "\t" final_species_taxid
}
EOF


        # Execute species lookup
        echo "  Running species name lookups..."
        awk -f species_lookup.awk "$NAMES_DMP" unique_species_names.txt > species_taxid_mapping.txt
        
        # Validate mapping file was created
        if [ ! -f "species_taxid_mapping.txt" ]; then
            echo "Error: Species mapping file not created"
            cd ..
            rm -rf "$TEMP_DIR_PART2"
            exit 1
        fi
        
        SPECIES_LOOKUP_SUCCESS=$(grep -v $'\tN/A$' species_taxid_mapping.txt | wc -l)
        SPECIES_SUCCESS_RATE=$(echo "scale=1; $SPECIES_LOOKUP_SUCCESS * 100 / $UNIQUE_SPECIES_COUNT" | bc -l 2>/dev/null || echo "unknown")
        
        echo "  Species lookup results: $SPECIES_LOOKUP_SUCCESS/$UNIQUE_SPECIES_COUNT found ($SPECIES_SUCCESS_RATE%)"
        
        # Step 4: Update the TSV with species-based taxids
        echo "Part II Step 4: Updating N/A rows with species-based taxids..."
        
        # Create AWK script to update N/A rows
        cat > update_na_rows.awk << 'EOF'
BEGIN {
    FS = OFS = "\t"
    
    # Load species name -> taxid mapping
    while ((getline line < "species_taxid_mapping.txt") > 0) {
        split(line, parts, "\t")
        if (length(parts) >= 2) {
            species_map[parts[1]] = parts[2]
        }
    }
    close("species_taxid_mapping.txt")
}

{
    # If this row has N/A in column 8 and valid species name in column 108
    if (NF >= 108 && $8 == "N/A" && $108 != "" && $108 != "N/A") {
        species_name = $108
        original_taxname = $7
        new_taxid = (species_name in species_map) ? species_map[species_name] : "N/A"
        
        if (new_taxid != "N/A") {
            print "DEBUG: Found Part II match for " original_taxname " (species name " species_name ")" > "part2_matches.log"
            $8 = new_taxid  # Update the taxid column
        } else {
            print "DEBUG: No Part II match for " original_taxname " (species name " species_name ")" > "part2_matches.log"    
        }
    }
    
    print $0
}
EOF

        # Apply the update and create final file
        echo "  Creating final output file..."
        zcat "$OUTPUT_FILE" | awk -f update_na_rows.awk | gzip > "$PART2_OUTPUT_FILE"
        
        # Validate final output was created
        if [ ! -f "$PART2_OUTPUT_FILE" ]; then
            echo "Error: Final output file not created"
            cd ..
            rm -rf "$TEMP_DIR_PART2"
            exit 1
        fi
        
        # Generate Part II statistics
        FINAL_NA_ROWS=$(zcat "$PART2_OUTPUT_FILE" | tail -n +2 | cut -f8 | grep -c '^N/A$')
        PART2_RECOVERED=$((NA_ROWS - FINAL_NA_ROWS))
        
        echo "Part II Results:"
        echo "  Original N/A rows: $NA_ROWS"
        echo "  Remaining N/A rows: $FINAL_NA_ROWS"
        echo "  Recovered by Part II: $PART2_RECOVERED"
        if [ "$NA_ROWS" -gt 0 ]; then
            PART2_SUCCESS_RATE=$(echo "scale=1; $PART2_RECOVERED * 100 / $NA_ROWS" | bc -l 2>/dev/null || echo "unknown")
            echo "  Part II recovery rate: $PART2_SUCCESS_RATE%"
        fi
        
        # Save Part II matches log before cleanup
        if [ -f "part2_matches.log" ]; then
            PART2_MATCHES_COUNT=$(wc -l < part2_matches.log)
            cp "part2_matches.log" "../part2_matches.log"
            echo "  Part II matches logged: $PART2_MATCHES_COUNT cases saved to part2_matches.log"
        fi
        
        cd ..
        rm -rf "$TEMP_DIR_PART2"
        ## combine part I & part II matches.log
        cat multiple_matches_and_partI_assignments.log part2_matches.log > log_multiple_matches_and_assignments.txt
        rm multiple_matches_and_partI_assignments.log
        rm part2_matches.log
        mv log_multiple_matches_and_assignments.txt log_multiple_matches_and_assignments_viruses.txt

        echo ""
        echo "Part II completed successfully!"
        echo "Final output file: $PART2_OUTPUT_FILE"
        
        # Update overall statistics
        OVERALL_FINAL_NA=$(zcat "$PART2_OUTPUT_FILE" | tail -n +2 | cut -f8 | grep -c '^N/A$')
        OVERALL_TOTAL=$(zcat "$PART2_OUTPUT_FILE" | tail -n +2 | wc -l)
        OVERALL_FOUND=$((OVERALL_TOTAL - OVERALL_FINAL_NA))
        OVERALL_SUCCESS_RATE=$(echo "scale=1; $OVERALL_FOUND * 100 / $OVERALL_TOTAL" | bc -l 2>/dev/null || echo "unknown")
        
        echo ""
        echo "Overall Final Statistics:"
        echo "  Total rows: $OVERALL_TOTAL"
        echo "  Final taxids found: $OVERALL_FOUND"
        echo "  Final N/A remaining: $OVERALL_FINAL_NA"
        echo "  Overall success rate: $OVERALL_SUCCESS_RATE%"
    fi
fi