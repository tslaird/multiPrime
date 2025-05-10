#!/bin/bash

# List of directories to process
#DIRS="dir1 dir2 dir3"
# Or use a glob pattern like this:
DIRS="input/EEE input/EHV input/Neospora_hughesi input/Sarcocystis_neurona input/Theileria_equi"

# Check if DIRS is empty
if [ -z "$DIRS" ]; then
  echo "No directories found"
  exit 1
fi

for dir in $DIRS; do
  if [ -d "$dir" ]; then
    echo "Processing directory: $dir"
    FULL_DIR=$(realpath "$dir")
    cd "$dir"
    pwd
    ALL_FA=$(ls fasta/all*.fa 2>/dev/null | head -n 1)
    ALL_FA_BASE=${ALL_FA##*/}
    echo $ALL_FA_BASE
    echo $ALL_FA
    YAML_FILE=$(ls *yaml 2>/dev/null | head -n 1)
    #echo "ALL_FA: $ALL_FA"
    #echo "YAML_FILE: $YAML_FILE"
    # Check if required files exist
    if [ -f "$ALL_FA" ] && [ -f "$YAML_FILE" ]; then
      SENTINEL_FILE="${ALL_FA}.processed"
      if [ ! -f "$SENTINEL_FILE" ] || [ "$ALL_FA" -nt "$SENTINEL_FILE" ]; then
        # convert fasta headers and seq to unique hash
        ../../scripts/hash_fasta.py -i $ALL_FA -o ${ALL_FA_BASE}.hash.fa -m ${ALL_FA_BASE}.hashmap
        # edit the config yaml
        #sed -i "/virus:\s*-/s/.*/  - ${ALL_FA%.fa}.hash/" $YAML_FILE
        # run snakemake multiPrime workflow
        #TODO figure out how to override and set config file values
        snakemake --configfile $YAML_FILE -s ../../multiPrime.py --cores 4 --resources mem_mb=60000 -C input_dir=$FULL_DIR results_dir=$FULL_DIR/results log_dir=$FULL_DIR/logs virus=${ALL_FA_BASE}.hash && touch .snakemake.done
        if [ -f "results/Core_primers_set/BWT_coverage/core_final_maxprimers_set.out" ] && [ -f ".snakemake.done" ]; then
            # extract amplicons
            ../../scripts/extractAmplicon.py -t results/Core_primers_set/BWT_coverage/core_final_maxprimers_set.out -f ${ALL_FA_BASE}.hash.fa -o ${ALL_FA_BASE}_ampicons_hash.fasta
            # map back to original headers
            ../../scripts/map_hash_fasta.py -m ${ALL_FA_BASE}.hashmap -i ${ALL_FA_BASE}_ampicons_hash.fasta -o ${ALL_FA_BASE}_ampicons.fasta
            # make final output table
            ../../scripts/build_table.py -t results/Core_primers_set/BWT_coverage/core_final_maxprimers_set.out -p results/Core_primers_set/core_final_maxprimers_set.fa -g results/Total_fa/${ALL_FA_BASE}.hash.format.fa -o ${ALL_FA_BASE}_primer_amplicon_table.csv -m ${ALL_FA_BASE}.hashmap
            # Update sentinel file
            touch "$SENTINEL_FILE"
        fi
      else
        echo "ALL_FA has not been updated, skipping processing of $dir"
      fi
    else
      echo "Required files not found in $dir, skipping"
    fi

    cd - > /dev/null  # Return to the previous directory
  else
    echo "Skipping $dir as it's not a directory"
  fi
done

echo "All directories processed"
