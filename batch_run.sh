#!/bin/bash

# List of directories to process
#DIRS="dir1 dir2 dir3"
# Or use a glob pattern like this:
DIRS="input/EEE"

# Check if DIRS is empty
if [ -z "$DIRS" ]; then
    echo "No directories found"
    exit 1
fi

for dir in $DIRS; do
    if [ -d "$dir" ]; then
        echo "Processing directory: $dir"
        cd "$dir"
        pwd
        FULL_DIR=$(realpath "$dir")
	ALL_FA=$(ls all*.fa 2>/dev/null | head -n 1)
	YAML_FILE=$(ls *yaml 2>/dev/null | head -n 1)
	echo "ALL_FA: $ALL_FA"
	echo "YAML_FILE: $YAML_FILE"
        # Check if required files exist
        if [ -f "$ALL_FA" ] && [ -f "$YAML_FILE" ]; then
            #convert fasta headers and seq to unique hash
            ../../scripts/hash_fasta.py -i $ALL_FA -o ${ALL_FA%.fa}.hash.fa -m ${ALL_FA%.fa}.hashmap
            # edit the config yaml
            #sed -i "/virus:\s*-/s/.*/  - ${ALL_FA%.fa}.hash/" $YAML_FILE
            # run snakemake multiPrime workflow
            #TODO figure out how to override and set config file values 
            snakemake --configfile $YAML_FILE -s ../../multiPrime.py --cores 4 --resources mem_mb=60000
            # extract amplicons
            ../../scripts/extractAmplicon.py -t results/Core_primers_set/BWT_coverage/core_final_maxprimers_set.out -f ${ALL_FA%.fa}.hash.fa -o ${ALL_FA%.fa}_ampicons_hash.fasta
            # map back to original headers
            ../scripts/map_hash_fasta.py -m ${ALL_FA%.fa}.hashmap -i ${ALL_FA%.fa}_ampicons_hash.fasta -o {ALL_FA%.fa}_ampicons.fasta
            # make final output table
            ../scripts/build_table.py -t results/Core_primers_set/BWT_coverage/core_final_maxprimers_set.out -p results/Core_primers_set/core_final_maxprimers_set.fa -g results/Total_fa/${ALL_FA%.fa}.hash.format.fa -o ${ALL_FA%.fa}_primer_amplicon_table.csv -m ${ALL_FA%.fa}.hashmap
        else
            echo "Required files not found in $dir, skipping"
        fi

        cd - > /dev/null  # Return to the previous directory
    else
        echo "Skipping $dir as it's not a directory"
    fi
done

echo "All directories processed"
