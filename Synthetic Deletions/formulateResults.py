import os
import sys
import pandas as pd

# Store the data from the VCF file in a pair containing the SVLEN and SVTYPE
def get_info(filename):
    data = []
    with open(filename, "r") as file:
        for line in file:
            if line[0] == "#":
                continue
            line = line.split("\t")
            info = line[7].split(";")
            svlen = 0
            svstart = int(line[1])
            svend = 0
            svtype = ""
            for i in info:
                if i[:5] == "SVLEN":
                    svlen = abs(int(i[6:]))
                if i[:6] == "SVTYPE":
                    svtype = i[7:]
                if i[:3] == "END":
                    svend = int(i[4:])
            data.append((svstart, svend, svlen, svtype))
    return data

# Function to count variants in each bin for each coverage
def count_variants_in_bins(vcf_data, bins):
    counts = {bin_range: 0 for bin_range in bins}
    for index, row in vcf_data.iterrows():
        svlen = abs(row['SVLEN'])  # Use absolute value of SVLEN
        for bin_range in bins:
            if bin_range[0] <= svlen <= bin_range[1]:
                counts[bin_range] += 1
                break  # Move to the next variant
    return counts

def tabulate_data():
    filecount = 0
    filesNotFound = []

    Data_Type = ["CCS"]
    zygousity = ["Homozygous", "Heterozygous"]
    Variant_Callers = ["SVIM", "Sniffles", "CuteSV"]
    Coverage = ["5x", "10x", "20x", "30x", "40x", "50x"]

    # Define the bins with ranges
    bins = [(50, 100), (300, 500), (501, 1000), (1500, 2000),
            (4000, 5000), (8000, 10000), (25000, 30000)]

    for zyg in zygousity:
        for data_type in Data_Type:
            for caller in Variant_Callers:
                summary_filename = f"{zyg}/Results/{caller}/{caller}_summary.csv"

                # Initialize DataFrame to store bin counts
                bin_counts_df = pd.DataFrame(columns=['Bin'] + Coverage)

                for coverage in Coverage:
                    counts_per_coverage = {bin_range: 0 for bin_range in bins}
                    for k in range(1,11):
                        filename = f"{zyg}/VCF_Files/{coverage}/modRand{coverage}{data_type}SVs{caller}_{k}.vcf"
                        if not os.path.exists(filename):
                            filesNotFound.append(filename)
                            filecount += 1
                            continue
                        else:
                            data = get_info(filename)
                            data.sort(key=lambda x: x[2])  # Sort by SV length
                            vcf_data = pd.DataFrame(data, columns=["Start", "End", "SVLEN", "SVTYPE"])
                            counts = count_variants_in_bins(vcf_data, bins)
                            for bin_range in bins:
                                counts_per_coverage[bin_range] += counts[bin_range]

                    for bin_range in bins:
                        bin_str = f"{bin_range[0]}-{bin_range[1]}"
                        row_exists = bin_counts_df['Bin'].eq(bin_str).any()
                        if row_exists:
                            bin_counts_df.loc[bin_counts_df['Bin'] == bin_str, coverage] = counts_per_coverage[bin_range]
                        else:
                            row = {'Bin': bin_str, coverage: counts_per_coverage[bin_range]}
                            bin_counts_df = pd.concat([bin_counts_df, pd.DataFrame([row])], ignore_index=True)

                os.makedirs(os.path.dirname(summary_filename), exist_ok=True)
                bin_counts_df.to_csv(summary_filename, index=False)

    if filecount == 0:
        print("All files found")
    else:
        print(f"Number of files not found: {filecount}")
        print("Files not found:")
        for file in filesNotFound:
            print(file)

tabulate_data()