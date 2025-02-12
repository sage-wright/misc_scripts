from Bio import SeqIO
import pysam
import warnings
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import groupby, count

warnings.filterwarnings("ignore")

def parse_gene_positions(genbank_file):
    """This function extracts all gene positions from a GenBank file.

    Args:
        genbank_file (File): the GBFF file from GenBank to parse

    Returns:
        Dictionary: the name of the gene as the key mapped to a value of [start, end] positions
    """
    record = SeqIO.read(genbank_file, "genbank")
    ref_gene_positions = {}

    for feature in record.features:
        if feature.type == "gene":
            start = int(feature.location.start)
            end = int(feature.location.end)
            gene_name = feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else 'unknown'
            ref_gene_positions[gene_name] = {"start": start, "end": end}

    return ref_gene_positions

def parse_primer_ranges(primer_file):
    """This function extracts the start/end positions from a BED-like file.

    Args:
        primer_file (File): A three-column tab-delimited file with the following columns: 0 - primer name, 1 - start position, 2 - end position

    Returns:
        Dictionary: the name of the primer as the key mapped to a value of [start, end] positions
    """ 
    primer_positions = {}

    with open(primer_file, "r") as infile:
        for line in infile:
            line = line.strip().split('\t')
            start = int(line[1])
            end = int(line[2])
            gene_name = line[0]
            primer_positions[gene_name] = {"start": start, "end": end}
    
    return primer_positions

def print_high_coverage_ranges(bam_file, target_positions):
    """This function determines how much of the range as indicated by the target_positions dictionary is covered by at least 20x coverage.

    Args:
        bam_file (File): The BAM file to determine coverage from
        target_positions (Dictionary): A dictionary with the following format: {name: {"start": start, "end": end}}
    """

    samfile = pysam.AlignmentFile(bam_file, "rb")

    for gene, positions in target_positions.items():
        start = positions["start"]
        end = positions["end"]
        high_coverage_ranges = []
        current_range = None

        for pileupcolumn in samfile.pileup("Chromosome", start, end): 
            if pileupcolumn.n >= 20:  # Check if coverage is higher than 20x
                if current_range is None:
                    current_range = [pileupcolumn.pos, pileupcolumn.pos]
                else:
                    current_range[1] = pileupcolumn.pos
            elif current_range is not None:
                high_coverage_ranges.append(tuple(current_range))
                current_range = None

        if current_range is not None:
            high_coverage_ranges.append(tuple(current_range))

        print(f"Gene: {gene}, Gene range: {start}-{end}, High coverage ranges: {high_coverage_ranges}")

    samfile.close()

def as_range(iterable):
    """_summary_This function takes an iterable and transforms it into a string range with the format start-stop.

    Args:
        iterable (List): A list of positions to transform into a range in sequential order (e.g., 1,2,3,4)

    Returns:
        String: The range coverted into a range representation (e.g., 1-4)
    """
    l = list(iterable)
    if len(l) > 1:
        return '{0}-{1}'.format(l[0], l[-1])
    else:
        return '{0}'.format(l[0])

def compute_common_high_coverage_ranges(bam_files, target_positions):
    """This function determines the common high coverage (>20x) ranges for a set of BAM files.

    Args:
        bam_files (List[File]): A list of BAM files to parse
        target_positions (Dictionary): A dictionary with the format of {name: {"start": start, "end": end}} where (start,end) is a tuple

    Returns:
        Dictionary: A dictionary with the format of {name: [(start, end),(start2,end2),...]}
    """
    common_high_coverage_ranges = {}
    ranges_dict = {}

    for bam_file in bam_files:
        samfile = pysam.AlignmentFile(bam_file, "rb")
        ranges_dict[bam_file] = {}

        for gene, positions in target_positions.items():
            start = positions["start"]
            end = positions["end"]
            high_coverage_ranges = []
            tuple_ranges = []
            current_range = None

            for pileupcolumn in samfile.pileup("Chromosome", start, end):  # Replace "H37Rv" with your reference sequence name
                if pileupcolumn.n >= 20:  # Check if coverage is higher than 20x
                    if current_range is None:
                        current_range_start = pileupcolumn.pos
                        current_range_end = pileupcolumn.pos
                        current_range = [current_range_start, current_range_end]
                    else:
                        current_range_end = pileupcolumn.pos
                        current_range[1] = current_range_end
                elif current_range is not None:
                    high_coverage_ranges.extend(range(current_range_start, current_range_end))
                    tuple_ranges.append(tuple([current_range_start, current_range_end]))
                    current_range = None

            if current_range is not None:
                high_coverage_ranges.extend(range(current_range_start, current_range_end))
                tuple_ranges.append(tuple([current_range_start, current_range_end]))

            ranges_dict[bam_file][gene] = tuple_ranges
            
            if gene not in common_high_coverage_ranges:
                common_high_coverage_ranges[gene] = set(high_coverage_ranges)
            else:
                common_high_coverage_ranges[gene] &= set(high_coverage_ranges)

        samfile.close()

    print("Common high coverage ranges:")
    
    #print(common_high_coverage_ranges)
    for gene, ranges in common_high_coverage_ranges.items():
        print(f"Gene: {gene}")
        print(','.join(as_range(g) for _, g in groupby(sorted(list(ranges)), lambda n, c=count(): n-next(c))))
        
    return ranges_dict

def ranges_to_dataframe(ranges_dict):
    """This function converts a dictionary of ranges into a pandas DataFrame.

    Args:
        ranges_dict (Dictionary): A dictionary with the format of {name: [(start, end),(start2,end2),...]} where (start,end) is a tuple

    Returns:
        pd.DataFrame: A DataFrame with the following columns: "BAM File", "Gene", "Position"
    """
    data = []

    for bam_file, gene_ranges in ranges_dict.items():
        for gene, ranges in gene_ranges.items():
            for coverage_range in ranges:
                for position in np.arange(coverage_range[0], coverage_range[1] + 1):
                    data.append({"BAM File": bam_file, "Gene": gene, "Position": position})

    return pd.DataFrame(data)

def dictionary_to_dataframe(target_positions):
    """This function converts a dictionary to DataFrame

    Args:
        target_positions (Dictionary): A dictionary with the following format: {name: {"start": start, "end": end}}

    Returns:
        pd.DataFrame: A DataFrame with three columns: "BAM File", "Gene", "Position"
    """
    df = pd.DataFrame()
    for gene, positions in target_positions.items():
        for position in np.arange(positions["start"], positions["end"] + 1):
            new_row = pd.DataFrame({"BAM File": ["H37Rv"], "Gene": [gene], "Position": [position]})
            df = pd.concat([df, new_row], ignore_index=True)

    return df

def plot_ranges_seaborn(df, df_ref, df2, bam_files):
    """This function plots two pd.DataFrames using seaborn, coloring the points based on the dataset they belong to.

    Args:
        df (pd.DataFrame): A DataFrame with three columns: "BAM File", "Gene", "Position"
        df2 (pd.DataFrame): A DataFrame with three columns: "BAM File", "Gene", "Position"
    """
    genes = df["Gene"].unique()

    for gene in genes:
        plt.figure(figsize=(10, 10))

        # Filter the dataframe for the current gene
        df_gene = df[df["Gene"] == gene]
        df2_gene = df2[df2["Gene"] == gene]
        df_ref_gene = df_ref[df_ref["Gene"] == gene]

        df_gene['dataset'] = "dataset1"
        df2_gene['dataset'] = "dataset2"
        df_ref_gene['dataset'] = "reference"

        combined = pd.concat([df_gene, df2_gene, df_ref_gene])
        # Plot the positions
        sns.scatterplot(data=combined, x="Position", y="BAM File", s=100, hue="dataset", palette=["black", "red", "green"], legend=False, edgecolor=None)

        #sns.scatterplot(data=df2_gene, x="Position", y="BAM File", s=100, color="red", legend=False, edgecolor=None)
        
        # Mark the positions in the plot where there's dots on all the samples
        for position in df_gene["Position"].unique():
            # ignore the reference

            #print([bam_file for bam_file in df_gene[df_gene["Position"] == position]["BAM File"].unique() if bam_file != "H37Rv"])
            gene_list = [bam_file for bam_file in df_gene[df_gene["Position"] == position]["BAM File"].unique() if bam_file != "H37Rv"]
            if sorted(gene_list) == sorted(bam_files):
                plt.axvline(position, color="grey", linestyle="-", linewidth=1, alpha=0.4)

        plt.xlabel('Position')
        plt.ylabel('BAM File')
        plt.title(f'High Coverage Ranges for {gene}')
        
        # Save the plot
        plt.savefig(f"high_coverage_ranges_{gene}.png")

def find_adjacent_ranges(lst):
    """This function determines if any ranges in a list are adjacent to each other and merges them.

    Args:
        lst (List): A list of positions to check for adjacency

    Returns:
        List: A list of lists where the first value is the start of the range and the second is the end of the range
    """
    if not lst:
        return []

    lst.sort()
    ranges = [[lst[0], lst[0]]]

    for num in lst[1:]:
        if num == ranges[-1][1] + 1:
            ranges[-1][1] = num
        else:
            ranges.append([num, num])

    return ranges

def get_common_positions(df, bam_files, target_positions):
    """This function determines the range of high coverage positions across all BAM files

    Args:
        df (pd.DataFrame): a DataFrame with three columns: "BAM File", "Gene", "Position"
        bam_files (List[File]): A list of BAM files to compare and process
        target_positions (Dictionary): A dictionary with the following format: {name: {"start": start, "end": end}}

    Returns:
        Dictionary: A dictionary with the following format: {name: [[start, end],[start2,end2],...]}
    """
    common_positions = {}

    for gene, positions in target_positions.items():
        common_positions[gene] = []

        for position in np.arange(positions["start"], positions["end"] + 1):
            
            gene_list = [bam_file for bam_file in df[df["Position"] == position]["BAM File"].unique() if bam_file != "H37Rv"]
            if sorted(gene_list) == sorted(bam_files):
                common_positions[gene].append(position)
    
    # merge the positions to the common_positions dictionary if they are adjacent to each other by updating max and min values
    for gene, positions in common_positions.items():
        common_positions[gene] = find_adjacent_ranges(positions)

    return common_positions


def __main__():
    # a list of all input BAM files (index file should be in same directory)
    bam_files = ["/home/sage_wright/github/misc_scripts/coverage_locations/bams/24TBEO-04-NY2P.bam", 
                "/home/sage_wright/github/misc_scripts/coverage_locations/bams/24TBEO-05-NY2P.bam",
                "/home/sage_wright/github/misc_scripts/coverage_locations/bams/24TBEO-15-NY2P.bam",
                "/home/sage_wright/github/misc_scripts/coverage_locations/bams/24TBEO-16-NY2P.bam",
                "/home/sage_wright/github/misc_scripts/coverage_locations/bams/24TBEO-17-NY2P.bam",
                "/home/sage_wright/github/misc_scripts/coverage_locations/bams/24TBEO-18-NY2P.bam",
                "/home/sage_wright/github/misc_scripts/coverage_locations/bams/24TBEO-26-NY2P.bam",
                "/home/sage_wright/github/misc_scripts/coverage_locations/bams/24TBEO-27-NY2P.bam"             
                ]
    # this file takes the format of primer_name\tstart\tend
    primer_range_file = "/home/sage_wright/github/misc_scripts/coverage_locations/primer-ranges-nymain.txt"
    # this file has the header of "BAM File"\t"Position"\t"Gene" and can be generated by the script `get_full_primer_overlap_range.py`
    # the input to that python script is the formatted*.tsv file; which is a formatted verion of primer-ranges*.tsv with primer_name\tstart\tend\ttarget_gene
    full_primer_range_file = "/home/sage_wright/github/misc_scripts/coverage_locations/full-primer-range-nymain.tsv"
    # the GenBank file for the reference genome
    genbank_file = "/home/sage_wright/github/misc_scripts/coverage_locations/H37Rv.gb"  

    # parse the GenBank file for gene positions
    ref_gene_positions = parse_gene_positions(genbank_file)

    # filter for just the genes of interest in the tNGS scheme
    # this is a list for the **full tNGS scheme**
    genes_of_interest = ["ahpC", "atpE", "eis", "embB", "embC", "embA", "ethA", 
                        "fabG1", "gyrA", "gyrB","inhA", "katG", "oxyR'", "pncA",
                        "rplC", "rpoB", "rpsL", "rrl", "rrs"] 
    genes_of_interest_nymain = ["ahpC", "eis", "embA", "embB", "embC", "ethA", 
                                "fabG1", "gyrA", "gyrB", "katG", "inhA", "oxyR'", 
                                "pncA", "rpoB", "rpsL", "rrs"]
    # "IS61110", "Rv0678/mmpR5" are not found in the GenBank file

    # make sure this is swapped when doing the main analysis
    ref_gene_positions = {gene: ref_gene_positions[gene] for gene in genes_of_interest_nymain}
    ref_gene_positions["Rv0678"] = {"start": 778990, "end": 779487}

    # map against the primer sequences
    primer_positions = parse_primer_ranges(primer_range_file)

    #ranges_dict = compute_common_high_coverage_ranges(bam_files, ref_gene_positions)
    ranges_dict = compute_common_high_coverage_ranges(bam_files, primer_positions)

    df = ranges_to_dataframe(ranges_dict)
    df_ref = dictionary_to_dataframe(ref_gene_positions)

    primers = pd.read_csv(full_primer_range_file, sep="\t", header=0)
    plot_ranges_seaborn(df, df_ref, primers, bam_files)

    common_positions = get_common_positions(df, bam_files, primer_positions)
    print(common_positions)

    for file in bam_files:
        print(f"\n{os.path.basename(file)}")
        print_high_coverage_ranges(file, primer_positions)
   
if __name__ == "__main__":
    __main__()