import pandas as pd
import argparse
import os
import sys

def get_input_args():
  parser = argparse.ArgumentParser(
    description = "This script anticipates concatenated laboratorian reports from tbp-parser where the first column is multiple repeated sample IDs that have unique suffixes. Limit of detection studies typically anticipate the same sample appearing multiple times, typically at different levels of depth or coverage. In order to determine at what coverage level the limit of detection is reached, this script will compare all samples with the same prefix and indicate at the suffix where LOD was obtained. Sample IDs should take the format of `name<suffix>`, where suffix can be any string but must end in an underscore and number that indicates the expected coverage level for the sample (sample01_1 or sample_dr_loci_20; coverage levels 1 and 20 respectively). The required input 'suffix file' lists each of the anticipated suffixes, one per line. The sample that is considered the 'truth' will need to be specified in the truth input variable. For more information, please reach out to support@theiagen.com or sage.wright@theiagen.com"
  )
  parser.add_argument('laboratorian_reports', type=str, help='path to the concatenated laboratorian reports')
  parser.add_argument('-o', '--output', type=str, help='prefix for the output sorted csv')
  parser.add_argument('suffix_file', type=str, help='path to the suffix file; one line per suffix')
  parser.add_argument('truth', type=str, help='the suffix that is considered the ground truth')
  parser.add_argument('-mcl', '--min_cov_lvl', type=int, required=False, help='the minimum down-sampling level to consider')
  parser.add_argument('-l', '--lims_only', required=False, action='store_true', help='only include genes for LIMS')
  return parser.parse_args()

def get_highest_fn_stats(header, results, stat_type):
  #current header for final mutation results file
  #results = the dictionary where (key = mutation) and (value = all stats for that mutation) for a single sample.
  #Find the highest down-sampling level/depth/read support/frequency where a mutation was found before a FN result.

  fn_stats = None
  stat_index = {}
  valid_stat_types = ["Down-sampling Level", "Sample Depth", "Sample Read Support", "Sample Frequency"]
  if stat_type not in valid_stat_types:
    raise ValueError(f"Invalid stat_type: {stat_type}. Must be one of {valid_stat_types}")

  for i, e in enumerate(header):
    # minus 1 because the first column is the sample name
    stat_index[e] = i - 1

  try:
    # Find the highest stat for the given stat_type
    fn_stats = max(
      (_ for _ in results.items() if (_[1][10] == "FN" and _[1][0] != "-")), key=lambda item: item[1][stat_index[stat_type]], default=None
    )
  except:
    pass

  if fn_stats:
    # sometimes there are more than one mutation with the same highest 'stat' value - report the first and the count of duplicates
    dupe_max_count = sum(1 for _ in results.values() if _[stat_index[stat_type]] == fn_stats[1][stat_index[stat_type]] and _[10] == "FN" and _[0] != "-") - 1
    # metrics = "{stat_type}", "Gene", ("AA-Mutation" or "NT-Mutation")
    metrics = (
      f"{fn_stats[1][stat_index[stat_type]]} "
      f"[{fn_stats[1][stat_index["Gene"]]} ~"
      f"{fn_stats[1][stat_index["AA-Mutation"]] if fn_stats[1][stat_index["AA-Mutation"]] != 'NA' else fn_stats[1][stat_index["NT-Mutation"]]}] "
      f"{"[+" + str(dupe_max_count) + " more]" if dupe_max_count >= 1 else ''}"
    )
  else:
    metrics = "-"

  summary_string = f"Highest {stat_type} where a mutation was found before FN result:\t{metrics} "
  return summary_string

def import_suffixes(suffix_file):
  with open(suffix_file, 'r') as file:
    suffixes = [suffix.rstrip() for suffix in file]
  return suffixes

def compare_against_truth(base_sample_name, suffixes, one_sample, truth_sample, min_cov_lvl=None):
  # extract all of the data for a single sample
  # for each mutation in the truth sample, check the other suffixes that found the mutation and keep track of when it is not found
  truth_sample = truth_sample.drop_duplicates(
    subset=[
      'sample_id',
      'tbprofiler_gene_name',
      'tbprofiler_variant_substitution_nt',
    ], ignore_index=True
  )
  one_sample = one_sample.drop_duplicates(
    subset=[
      'sample_id',
      'tbprofiler_gene_name',
      'tbprofiler_variant_substitution_nt',
    ], ignore_index=True
  )
  if min_cov_lvl:
    one_sample = one_sample[one_sample["sample_id"].str.split("_").str[-1].astype(int) >= min_cov_lvl]

  # Ordered by (descending) coverage levels.
  sorted_samples = []
  for suffix in suffixes:
    sorted_samples.append(f"{base_sample_name}{suffix}")
  sorted_samples = sorted(sorted_samples, key=lambda x: int(x.split('_')[-1]), reverse=True)

  # this dictionary has the mutation as the key and stores the highest coverage level before a mutation was last found (along with read support,frequency, etc.)
  mutation_dictionary = {}

  # All mutations found in the truth sample/reference
  truth_mutations = truth_sample[["tbprofiler_variant_substitution_nt", "read_support"]]
  truth_mutations = truth_mutations.drop_duplicates(subset=['tbprofiler_variant_substitution_nt'], ignore_index=True)

  for mutation in truth_mutations["tbprofiler_variant_substitution_nt"]:
    aa_mutation = truth_sample[truth_sample["tbprofiler_variant_substitution_nt"] == mutation]["tbprofiler_variant_substitution_aa"].values[0]
    gene_name = truth_sample[truth_sample["tbprofiler_variant_substitution_nt"] == mutation]["tbprofiler_gene_name"].values[0]
    truth_depth = int(float(truth_sample[truth_sample["tbprofiler_variant_substitution_nt"] == mutation]["depth"].values[0]))
    truth_freq = float(truth_sample[truth_sample["tbprofiler_variant_substitution_nt"] == mutation]["frequency"].values[0])
    truth_read_support = int(float(truth_mutations[truth_mutations["tbprofiler_variant_substitution_nt"] == mutation]["read_support"].values[0]))

    true_positives = []
    for sample_id in sorted_samples:
      # Checking if mutation exists in the subsamples dataframe (looking for True positives). If it does, append True, else False. see below.
      true_positives.append(True if mutation in one_sample[one_sample["sample_id"] == sample_id]["tbprofiler_variant_substitution_nt"].values else False)

    if all(true_positives):
      pass

    elif any(true_positives) and not all(true_positives):
      # The index of the last "False" value in this list will correspond to the highest coverage level where the mutation was NOT found.
      # Catches situations where the mutation disappears and then reappears - we want to count LOD at the coverage level right before it disappears for the first time.
      first_false_index = true_positives.index(False)

      # This is the subsample that has the lowest coverage level where the mutation was found (while considering disappearances).
      # Will definitely have a False value by this point - labeling the subsample as FN (even though its technically TP
      sample_choice = sorted_samples[first_false_index - 1]
      coverage_level = int(sample_choice.split("_")[-1])
      sample = one_sample[one_sample["sample_id"].str.endswith(sample_choice)]
      sample_depth = int(float(sample[sample["tbprofiler_variant_substitution_nt"] == mutation]["depth"].values[0]))
      sample_freq = float(sample[sample["tbprofiler_variant_substitution_nt"] == mutation]["frequency"].values[0])
      sample_read_support = int(float(sample[sample["tbprofiler_variant_substitution_nt"] == mutation]["read_support"].values[0]))
      mutation_dictionary[mutation] = [coverage_level, gene_name, mutation, aa_mutation, sample_depth, sample_read_support, sample_freq, truth_depth, truth_read_support, truth_freq, "FN"]

    elif not any(true_positives):
      mutation_dictionary[mutation] = ["-", gene_name, mutation, aa_mutation, "-", "-", "-", truth_depth, truth_read_support, truth_freq, "FN"]

    else:
      print("Something went wrong. You shouldn't be here.")
      breakpoint()

  # If the mutation is found in the subsample, but not in the reference, it is a false positive.
  false_positives = one_sample[~one_sample['tbprofiler_variant_substitution_nt'].isin(truth_mutations["tbprofiler_variant_substitution_nt"])]
  if not false_positives.empty:
    for i, fp_sample in false_positives.iterrows():
      fp_mutation = fp_sample["tbprofiler_variant_substitution_nt"]
      fp_aa_mutation = fp_sample["tbprofiler_variant_substitution_aa"]
      fp_gene_name = fp_sample["tbprofiler_gene_name"]
      fp_coverage_level = int(fp_sample["sample_id"].split("_")[-1])
      fp_depth = int(float(fp_sample["depth"]))
      fp_read_support = int(float(fp_sample["read_support"]))
      fp_freq = float(fp_sample["frequency"])

      # Add highest coverage FP subsample to the mutation_dictionary
      if fp_mutation not in mutation_dictionary:
        mutation_dictionary[fp_mutation] = [fp_coverage_level, fp_gene_name, fp_mutation, fp_aa_mutation, fp_depth, fp_read_support, fp_freq, "-", "-", "-", "FP"]
      elif fp_coverage_level > mutation_dictionary[fp_mutation][0]:
        mutation_dictionary[fp_mutation] = [fp_coverage_level, fp_gene_name, fp_mutation, fp_aa_mutation, fp_depth, fp_read_support, fp_freq, "-", "-", "-", "FP"]
  return mutation_dictionary

# def sort_df(original_df, prefix):
#   original_df["coverage_level"] = original_df["sample_id"].str.split("_").str[-1].replace("full", 999).apply(lambda x: int(x))
#   original_df["sample_name"] = original_df["sample_id"].str.split("_").str[0]
#   #print(original_df["coverage_level"].unique())
#   sorted_df = original_df.sort_values(by=["sample_name", "tbprofiler_variant_substitution_nt", "antimicrobial", "coverage_level"], ascending=[True, True, True, False], ignore_index=True)
#   filename = prefix + ".sorted_laboratorian_reports.csv"
#   sorted_df.to_csv(filename, index=False)

def main():
  options = get_input_args()
  # import the combined laboratorian reports
  original_df = pd.read_csv(options.laboratorian_reports, na_filter=False)

  # remove all WT mutations
  original_df = original_df[original_df["tbprofiler_variant_substitution_nt"] != "WT"]
  if options.lims_only:
    GENES_FOR_LIMS_WGS = [
      "atpE", "eis", "embA", "embB", "ethA", "fabG1", "gyrA",
      "gyrB", "inhA", "katG", "mmpL5", "mmpS5", "pepQ",
      "pncA", "rplC", "rpoB", "rrl", "rrs", "Rv0678", "tlyA"
    ]
    original_df = original_df[original_df["tbprofiler_gene_name"].isin(GENES_FOR_LIMS_WGS)]
  #sort_df(original_df, options.output)

  suffixes = import_suffixes(options.suffix_file)

  # get list of sample IDs without any and all suffixes
  all_sample_ids = original_df["sample_id"]
  for suffix in suffixes:
    all_sample_ids = all_sample_ids.str.removesuffix(suffix)
  all_sample_ids = all_sample_ids.str.removesuffix(options.truth)

  sample_ids = set(all_sample_ids.unique())
  #print("\n".join(sample_ids))
  # for each sample ID, extract the data and compare the suffixes
  lowest_coverage_levels = {}

  for sample in sample_ids:
    one_sample = original_df[original_df["sample_id"].str.startswith(sample)]
    # compare to the truth sample
    truth_sample = one_sample[one_sample["sample_id"].str.endswith(options.truth)]
    one_sample = one_sample[~one_sample["sample_id"].str.endswith(options.truth)]
    lowest_coverage_levels[sample] = compare_against_truth(sample, suffixes, one_sample, truth_sample, options.min_cov_lvl)

  # sort by sample name first.
  lowest_coverage_levels = dict(sorted(lowest_coverage_levels.items()))

  # sort all the mutation_dictionaries by coverage level
  for sample in lowest_coverage_levels:
    lowest_coverage_levels[sample] = dict(
      sorted(
        lowest_coverage_levels[sample].items(),
        key=lambda item: (item[1][0], item[1][1].lower()) if item[1][0] != "-" else (float("inf"), item[1][1].lower())
      )
    )

  # Loop through sorted dictionary and print out the results.
  with open(options.output, "w") as f:
    header = ["Sample", "Down-sampling Level", "Gene", "NT-Mutation", "AA-Mutation", "Sample Depth", "Sample Read Support", "Sample Frequency", "Reference Depth", "Reference Read Support", "Reference Frequency", "Status"]
    print(",".join(header), file=f)

    for sample, results in lowest_coverage_levels.items():
      print(f"{sample}", file=f)
      print(f"Sample: {sample}")

      fn_count = 0
      fp_count = 0
      for mutation, stats in results.items():
        if stats[10] == "FN":
          fn_count += 1
        elif stats[10] == "FP":
          fp_count += 1
        #print to file all FPs and FNs (mutation_results.csv)
        print(' ,' + ','.join(str(_) for _ in stats), file=f)
      print(f"Total FPs: {fp_count}")
      print(f"Total FNs: {fn_count}")

      print(get_highest_fn_stats(header, results, "Down-sampling Level").expandtabs(8))
      print(get_highest_fn_stats(header, results, "Sample Depth").expandtabs(20))
      print(get_highest_fn_stats(header, results, "Sample Read Support").expandtabs(8))
      print(get_highest_fn_stats(header, results, "Sample Frequency").expandtabs(20))

      print("\n", file=f)
      print("\n")

if __name__ == "__main__":
  main()
