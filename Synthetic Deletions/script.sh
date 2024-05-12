#!/bin/bash

set -e

# Redirect all output to output.txt
exec &> logs.txt

# coverage=("5x" "10x")
coverage=("5x" "10x" "20x" "30x" "40x" "50x")
# zygosity=("Heterozygous")
zygosity=("Homozygous" "Heterozygous")
VariantCaller=("Sniffles" "SVIM" "CuteSV")
DataType=("CCS")
# DataType=("CCS" "CLR" "ONT")
for data_type in "${DataType[@]}"; do
    directories+=("Raw_Seq_Files_$data_type")
    directories+=("Alignment_Files_$data_type")
    directories+=("Deletions_$data_type")
done

echo "Creating Directories"

# Create directories for each zygosity, directory type, and coverage level
for zyg in "${zygosity[@]}"; do
    mkdir -p "$zyg"
    mkdir -p "$zyg/Results"
    for data_type in "${DataType}";do
        for vc in "${VariantCaller[@]}"; do
            mkdir -p "$zyg/Results/$data_type/$vc"
        for dir in "${directories[@]}"; do
            mkdir -p "$zyg/$dir"
            mkdir -p "$zyg/$dir/$data_type"
            for cov in "${coverage[@]}"; do
                mkdir -p "$zyg/$dir/$cov"
                done
            done
        done
    done
done

echo "Directories created"
echo "Starting the pipeline"

NUM_ITERATIONS=10

for ((i=1; i<=$NUM_ITERATIONS; i++)); do
    echo "--------------------------------"
    echo "Iteration $i"

    python3 genModSeq.py
    echo "--------------------------------"

    for zyg in "${zygosity[@]}"; do
        for data_type in "${DataType}"; do
            echo "Starting preprocessing for $zyg" 
            cd "$zyg" || exit 1
            
            chmod +x ../PreProcessingCommands.sh
            echo "Starting preprocessing for $zyg"
            .././PreProcessingCommands.sh
            echo "Preprocessing completed for $zyg"
            echo "--------------------------------"
            
            cd Deletions_$data_type || exit 1
            chmod +x ../../SVCallerCommands.sh
            echo "Starting SV calling for $zyg"
            ../.././SVCallerCommands.sh $i
            echo "SV calling completed for $zyg"
            echo "--------------------------------"
            cd .. || exit 1
            cd .. || exit 1

            echo "Completed the $i-th iteration for $zyg $data_type data"
        done
    done

    echo "--------------------------------"
done

# remove all the content in the folders Alignment_Files from both Homozygous and Heterozygous
for zyg in "${zygosity[@]}"; do
    for data_type in "${DataType}"; do
        rm -r $zyg/Alignment_Files_$data_type/*
    done
done

echo "Formulating and combining the results"
python3 formulateResults.py
echo "Results organized"

echo "Pipeline completed"
