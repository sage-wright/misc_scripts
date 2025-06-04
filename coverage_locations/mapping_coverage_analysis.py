from Bio.SeqFeature import SeqFeature, FeatureLocation
from collections import defaultdict
import pysam
import os
import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Coverage assessment script")
    parser.add_argument("-b", "--bam_dir", required=True, help="Directory containing BAM files")
    parser.add_argument("-p", "--primer_range_file", required=True, help="Input ('BED-like') primer range file")
    parser.add_argument("-a", "--assay", required=True, choices=["cdc", "ny3d", "nymain", "cdc_full"], help="Assay type (must be one of: 'cdc', 'ny3d', 'nymain')")
    parser.add_argument("-q", "--min_base_quality", type=int, default=0, help="Bases below the minimum quality will not be counted (default: 0)")
    parser.add_argument("-m", "--min_mapping_quality", type=int, default=0, help="Only use reads above a minimum mapping quality. (default: 0)")
    parser.add_argument("--ignore_overlaps", action="store_true", help="Ignore overlapping primer regions (default: False)")
    
    return parser.parse_args()

def overlaps(loc1, loc2):
    """This function checks if two loci overlap.

    Args:
        loc1 (FeatureLocation): A FeatureLocation representing the start and end positions of the first locus
        loc2 (FeatureLocation): A FeatureLocation representing the start and end positions of the second locus

    Returns:
        bool: True if the loci overlap, False otherwise
    """
    overlap = (min(loc1.end, loc2.end) - max(loc1.start, loc2.start))
    return overlap > 0

def get_split_regions(loc1, loc2):
    """This function splits two loci into overlapping and non-overlapping regions.
    Args:
        loc1 (FeatureLocation): A FeatureLocation representing the start and end positions of the first locus
        loc2 (FeatureLocation): A FeatureLocation representing the start and end positions of the second locus
    Returns:
        loci1 (dict): A dictionary with keys "overlapping" and "non_overlapping" containing lists of FeatureLocations
        loci2 (dict): A dictionary with keys "overlapping" and "non_overlapping" containing lists of FeatureLocations
    """
    overlap_start = max(loc1.start, loc2.start)
    overlap_end = min(loc1.end, loc2.end)

    # If no overlap, return entire regions as non-overlapping
    if overlap_start >= overlap_end:
        return (
            {"overlapping": [], "non_overlapping": [loc1]},
            {"overlapping": [], "non_overlapping": [loc2]}
        )

    def split(loc):
        parts = {"overlapping": [], "non_overlapping": []}
        parts["overlapping"].append(FeatureLocation(overlap_start, overlap_end))
        if loc.start < overlap_start:
            parts["non_overlapping"].append(FeatureLocation(loc.start, overlap_start))
        if loc.end > overlap_end:
            parts["non_overlapping"].append(FeatureLocation(overlap_end, loc.end))
        return parts

    return split(loc1), split(loc2)

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

            feature = SeqFeature(FeatureLocation(start, end))
            primer_positions[gene_name] = feature

    # sort primer positions by start position
    primer_positions = dict(sorted(primer_positions.items(), key=lambda kv: int(kv[1].location.start)))
    return primer_positions

def parse_overlapping_regions(bam_file, primer_positions, min_base_quality, min_mapping_quality):
    """This function extracts the read names from non-overlapping regions of the primer positions.

    Args:
        bam_file (File): A BAM file - alignment file between H37Rv reference genome and the sample
        primer_positions (dict): A dictionary with the name of the primer as the key mapped to a value of [start, end] positions
    Returns:
        Dictionary: nested dictionary keyed by the name of the primer (gene) and mapped to a list of read names
        ex) primer_specific[gene]["read_names"] = list of read names in the primer region
    """

    # check if the BAM file is indexed
    bam = pysam.AlignmentFile(bam_file)
    chrom = bam.get_reference_name(0)
    if not bam.check_index():
        bam_index = pysam.IndexedReads(bam)

    primer_specific = defaultdict(lambda: defaultdict(list))

    for i, (gene, feature) in enumerate(primer_positions.items()):

        # check if any primers overlap each other and get those regions
        for j, (gene2, feature2) in enumerate(primer_positions.items()):
            if i >=j:
                continue
            if overlaps(feature.location, feature2.location):
                # get the overlapping regions
                loci, loci2 = get_split_regions(feature.location, feature2.location)

                # gene 1
                primer_specific[gene]["ovl_region"] += (loci["overlapping"])

                if primer_specific[gene]["non_ovl_region"]:
                    new_non_ovl_region = FeatureLocation(
                        max(loci["non_overlapping"][0].start, primer_specific[gene]["non_ovl_region"][0].start),
                        min(loci["non_overlapping"][0].end, primer_specific[gene]["non_ovl_region"][0].end),
                    )
                    primer_specific[gene]["non_ovl_region"] = [new_non_ovl_region]
                else:
                    primer_specific[gene]["non_ovl_region"] = loci["non_overlapping"]

                primer_specific[gene]["ovl_gene"].append(gene2)
                primer_specific[gene]["unique_reads"] = set()

                # gene 2
                primer_specific[gene2]["ovl_region"] += loci2["overlapping"]

                if primer_specific[gene2]["non_ovl_region"]:
                    new_non_ovl_region = SeqFeature(FeatureLocation(
                        max(loci2["non_overlapping"][0].start, primer_specific[gene2]["non_ovl_region"][0].start),
                        min(loci2["non_overlapping"][0].end, primer_specific[gene2]["non_ovl_region"][0].end),
                    ))
                    primer_specific[gene2]["non_ovl_region"] = [new_non_ovl_region]
                else:
                    primer_specific[gene2]["non_ovl_region"] = loci2["non_overlapping"]

                primer_specific[gene2]["ovl_gene"].append(gene)
                primer_specific[gene2]["unique_reads"] = set()

    # currently does not consider complex overlaps with more than two genes overlapping the same region.

    for gene, feature in primer_positions.items():
        # process genes with identified non-overlapping regions
        if gene in primer_specific and primer_specific[gene]["non_ovl_region"]:
            region = primer_specific[gene]["non_ovl_region"][0]
            non_ovl_start = region.start
            non_ovl_end = region.end

            for rec in bam.pileup(
                contig=chrom,
                start=int(non_ovl_start),
                end=int(non_ovl_end),
                min_base_quality=min_base_quality,
                min_mapping_quality=min_mapping_quality,
                truncate=True,
                until_eof=True,
            ):
                primer_specific[gene]["unique_reads"].update(rec.get_query_names())

    return primer_specific

def parse_bam(bam_file, primer_positions, primer_specific, min_base_quality, min_mapping_quality):
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
    chrom = bam.get_reference_name(0)
    if not bam.check_index():
        bam_index = pysam.IndexedReads(bam)
    
    # initialize stats dictionary
    for i, (gene, feature) in enumerate(primer_positions.items()):
        stats[gene] = {}
        stats[gene]["avg"] = 0
        stats[gene]["med"] = 0
        stats[gene]["min"] = 0
        stats[gene]["max"] = 0

    for gene, feature in primer_positions.items():
        primer_start = feature.location.start
        primer_end = feature.location.end
        primer_len = len(feature)

        # initialize min_read dictionary with specific coverage levels (reset for each primer/gene)
        for cov in cov_lvl:
            min_read[cov] = []

        primer_pos_list = []
        for rec in bam.pileup(
            contig=chrom,
            start=int(primer_start),
            end=int(primer_end),
            min_base_quality=min_base_quality,
            min_mapping_quality=min_mapping_quality,
            truncate=True,
            until_eof=True,
        ):
            if gene in primer_specific and primer_specific[gene]["unique_reads"]:
                valid_reads = [
                    read for read in rec.get_query_names()
                    if read in primer_specific[gene]["unique_reads"]
                ]
                num_reads = len(valid_reads)
            else:
                num_reads = rec.get_num_aligned()

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
    ignore_overlaps = args.ignore_overlaps
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
    elif assay == "cdc_full":
      genes_of_interest = [
        "ahpC", "atpE", "eis", "embB", "fabG", "gyrA", "gyrB",
        "inhA", "katG1", "katG2", "katG3", "katG4", "pncA",
        "rplC_1", "rplC_2", "rpoB170", "rpoB-RRDR",
        "rrl_1", "rrl_2", "rrs", "rv0678",
      ]
    else:
        raise ValueError("Assay type must be one of: 'cdc', 'ny3d', 'nymain'")

    primer_positions = parse_primer_ranges(primer_range_file)

    for bam in sorted(bam_files):
        if not ignore_overlaps:
            ovl_regions = parse_overlapping_regions(bam, primer_positions, min_base_quality, min_mapping_quality)
        else:
            ovl_regions = defaultdict(lambda: {"ovl_region": [], "non_ovl_region": [], "ovl_gene": [], "unique_reads": set()})
        samples[bam] = parse_bam(bam, primer_positions, ovl_regions, min_base_quality, min_mapping_quality)

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