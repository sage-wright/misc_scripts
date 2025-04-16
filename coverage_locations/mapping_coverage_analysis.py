import pysam
import os
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Coverage assessment script")
    parser.add_argument("-b", "--bam_dir", required=True, help="Directory containing BAM files")
    parser.add_argument("-p", "--primer_range_file", required=True, help="Input ('BED-like') primer range file")
    parser.add_argument("-a", "--assay", required=True, choices=["cdc", "ny3d", "nymain"], help="Assay type (must be one of: 'cdc', 'ny3d', 'nymain')")
    parser.add_argument("-q", "--min_base_quality", type=int, default=0, help="Bases below the minimum quality will not be counted (default: 0)")
    parser.add_argument("-m", "--min_mapping_quality", type=int, default=0, help="Only use reads above a minimum mapping quality. (default: 0)")
    
    return parser.parse_args()

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

def parse_bam(bam_file, primer_positions, min_base_quality, min_mapping_quality):
    """This function calculates the average, median, min, and max coverage for each primer region in a BAM file.

    Args:
        bam_file (File): A BAM file - alignment file between H37Rv reference genome and the sample
        primer_positions (dict): A dictionary with the name of the primer as the key mapped to a value of [start, end] positions
    Returns:
        Dictionary: nested dictionary keyed by the name of the primer (gene) and mapped to a count of reads based on a specific metric
        ex) stats[gene]["avg"] = average number of reads in the primer region
        ex) stats[gene]["50x"] = percentage of reads with coverage >= 50x
    """
    min_read = {}
    stats = {}
    cov_lvl = [50, 100, 200, 500, 1000]
    sample_name = os.path.basename(bam_file).split(".")[0]
    header = [
        "Loci", "Coordinates", "Average loci coverage", "Median loci coverage", "Minimum loci coverage", "Maximum loci coverage",
        "Query Coverage 50x", "Query Coverage 100x", "Query Coverage 200x", "Query Coverage 500x", "Query Coverage 1000x"
    ]
    print(f"{sample_name}")
    print(f"{'\t'.join(header)}")
    
    # check if the BAM file is indexed
    bam = pysam.AlignmentFile(bam_file)
    if not bam.check_index():
        bam_index = pysam.IndexedReads(bam)
    
    # sort primer positions by start position
    primer_positions = sorted(primer_positions.items(), key=lambda x: x[1]["start"])
    
    # initialize stats dictionary
    for i, coords in enumerate(primer_positions):
        gene = coords[0]
        stats[gene] = {}
        stats[gene]["avg"] = 0
        stats[gene]["med"] = 0
        stats[gene]["min"] = 0
        stats[gene]["max"] = 0
        
    for i, coords in enumerate(primer_positions):
        gene = coords[0]
        primer_start = primer_positions[i][1]["start"]
        primer_end = primer_positions[i][1]["end"]
        primer_len = abs(primer_end - primer_start) + 1

        # initialize min_read dictionary with specific coverage levels (reset for each primer/gene)
        for cov in cov_lvl:
            min_read[cov] = []

        primer_pos_list = []
        for rec in bam.pileup(min_base_quality=min_base_quality, min_mapping_quality=min_mapping_quality, until_eof=True):
            pos = rec.reference_pos
            num_reads = rec.get_num_aligned()
            if pos < primer_start:
                continue
            elif pos > primer_end:
                break
            elif pos >= primer_start and pos <= primer_end:
                primer_pos_list.append(num_reads)
                for cov in cov_lvl:
                    if num_reads >= cov:
                        min_read[cov].append(num_reads)

        if len(primer_pos_list) > 0:            
            stats[gene]["avg"] = np.mean(primer_pos_list)
            stats[gene]["med"] = np.median(primer_pos_list)
            stats[gene]["min"] = np.min(primer_pos_list)
            stats[gene]["max"] = np.max(primer_pos_list)
            
        stats[gene]["50x"] = (len(min_read[50])/primer_len) * 100
        stats[gene]["100x"] = (len(min_read[100])/primer_len) * 100
        stats[gene]["200x"] = (len(min_read[200])/primer_len) * 100
        stats[gene]["500x"] = (len(min_read[500])/primer_len) * 100
        stats[gene]["1000x"] = (len(min_read[1000])/primer_len) * 100

        deliverables = [
            gene, 
            f"[{primer_start},{primer_end}]", f"{stats[gene]['avg']:.2f}", f"{stats[gene]['med']}", f"{stats[gene]['min']}", f"{stats[gene]['max']}",
            f"{stats[gene]['50x']:.2f}", f"{stats[gene]['100x']:.2f}", f"{stats[gene]['200x']:.2f}", f"{stats[gene]['500x']:.2f}", f"{stats[gene]['1000x']:.2f}"
        ]
        print(f"{'\t'.join([str(x) for x in deliverables])}")

    # Calculate average coverage across all loci
    final_avg = np.mean([stats[gene]["avg"] for gene in stats])
    final_med = np.median([stats[gene]["avg"] for gene in stats])
    final_min = np.min([stats[gene]["avg"] for gene in stats])
    final_max = np.max([stats[gene]["avg"] for gene in stats])

    print(f"Average coverage across all loci: {final_avg:.2f}")
    print(f"Median among average loci coverages: {final_med:.2f}")
    print(f"Minimum among average loci coverages: {final_min:.2f}")
    print(f"Maximum among average loci coverages: {final_max:.2f}")
    
    nl = "\n"
    print(f"{nl}")
    
    return stats	
  
def __main__():
    args = parse_args()
    bam_files = [os.path.join(args.bam_dir, f) for f in os.listdir(args.bam_dir) if f.endswith(".bam")]
    primer_range_file = args.primer_range_file
    min_base_quality = args.min_base_quality
    min_mapping_quality = args.min_mapping_quality
    assay = args.assay
    
    samples = {}

    genes_of_interest = []
    if assay == "cdc":
        genes_of_interest = [
        	"gyrB", "gyrA", "rpoB-170", "rpoB-RRDR",
        	"inhA", "fabG1","katG", "pncA",
        ]
    elif assay == "ny3d":
        genes_of_interest = [
        	"rpoB-RD1", "rpoB-RD2", "katG-RD1", "katG-RD2",
        	"gyrA-RD1", "gyrA-RD2", "gyrB", "inhA",
        ]
    elif assay == "nymain":
        genes_of_interest = [
        	"gyrB", "gyrA", "rpoB", "rpsL", "rrs", "inhA", 
        	"fabG1","katG", "pncA", "eis", "ahpC", "oxyR'", 
        	"embC", "embA", "embB", "ethA",
        ]
    else:
        raise ValueError("Assay type must be one of: 'cdc', 'ny3d', 'nymain'")
     
    # genes_of_interest_nymain = [
    #     "gyrB", "gyrA", "rpoB", "rpsL", "rrs", "inhA", 
    #     "fabG1","katG", "pncA", "eis", "ahpC", "oxyR'", 
    #     "embC", "embA", "embB", "ethA",
    # ]
    # genes_of_interest_cdc = [
    #     "gyrB", "gyrA", "rpoB-170", "rpoB-RRDR",
    #     "inhA", "fabG1","katG", "pncA",  
    # ]
    # genes_of_interest_nys3drug = [
    #     "rpoB-RD1", "rpoB-RD2", "katG-RD1", "katG-RD2",
    #     "gyrA-RD1", "gyrA-RD2", "gyrB", "inhA"
    # ]

    primer_positions = parse_primer_ranges(primer_range_file)

    for bam in sorted(bam_files):
        samples[bam] = parse_bam(bam, primer_positions, min_base_quality, min_mapping_quality)

    # calculate average, median, min, and max coverage for each locus across all samples
    locus_coverage_stats = {}
    for gene in genes_of_interest:
        locus_coverage_stats[gene] = {
            "avg": [],
            "med": [],
            "min": [],
            "max": []
        }

    for bam, stats in samples.items():
        for gene, gene_stats in stats.items():
            locus_coverage_stats[gene]["avg"].append(gene_stats["avg"])
            locus_coverage_stats[gene]["med"].append(gene_stats["med"])
            locus_coverage_stats[gene]["min"].append(gene_stats["min"])
            locus_coverage_stats[gene]["max"].append(gene_stats["max"])

    print("Across all samples:")
    print("Gene\tAverage Coverage\tMedian Coverage\tMinimum Coverage\tMaximum Coverage")
    for gene, coverage_stats in locus_coverage_stats.items():
        avg_coverage = np.mean(coverage_stats["avg"])
        med_coverage = np.median(coverage_stats["med"])
        min_coverage = np.min(coverage_stats["min"])
        max_coverage = np.max(coverage_stats["max"])

        print(f"{gene}\t{avg_coverage:.2f}\t{med_coverage:.2f}\t{min_coverage:.2f}\t{max_coverage:.2f}")
   
if __name__ == "__main__":
    __main__()