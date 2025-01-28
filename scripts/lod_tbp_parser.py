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
  return parser.parse_args()

def import_suffixes(suffix_file):
  with open(suffix_file, 'r') as file:
    suffixes = [suffix.rstrip() for suffix in file]
  return suffixes

def compare_against_truth(one_sample, truth_sample):
  # extract all of the data for a single sample
  # for each mutation in the truth sample, check the other suffixes that found the mutation and keep track of when it is not found
  truth_sample = truth_sample.drop_duplicates(
    subset=[
      'sample_id',
      'tbprofiler_gene_name',
      'tbprofiler_variant_substitution_nt',
    ], ignore_index=True
  )
  truth_mutations = truth_sample[["tbprofiler_variant_substitution_nt", "read_support"]]
  truth_mutations = truth_mutations.drop_duplicates(subset=['tbprofiler_variant_substitution_nt'], ignore_index=True)
  # this dictionary has the mutation as the key and the LOWEST coverage level where the mutation was found and the mutation's read support at that coverage level
  mutation_dictionary = {}
  for mutation in truth_mutations["tbprofiler_variant_substitution_nt"]:
    gene_name = truth_sample[truth_sample["tbprofiler_variant_substitution_nt"] == mutation]["tbprofiler_gene_name"].values[0]
    truth_depth = int(truth_sample[truth_sample["tbprofiler_variant_substitution_nt"] == mutation]["depth"].values[0])
    truth_freq = float(truth_sample[truth_sample["tbprofiler_variant_substitution_nt"] == mutation]["frequency"].values[0])
    truth_read_support = int(truth_mutations[truth_mutations["tbprofiler_variant_substitution_nt"] == mutation]["read_support"].values[0])

    # Ordering by coverage levels.
    sorted_samples = sorted(one_sample["sample_id"].unique(), key=lambda x: int(x.split('_')[-1]))

    true_positives = []
    for suffix in sorted_samples:
      # Checking if mutation exists in the subsamples dataframe (looking for True positives). If it does, append True, else False. see below.
      true_positives.append(True if mutation in one_sample[one_sample["sample_id"].str.endswith(suffix)]["tbprofiler_variant_substitution_nt"].values else False)

      # If the mutation is found in the subsample, but not in the reference, it is a false positive.
      temp_sample = one_sample[one_sample["sample_id"].str.endswith(suffix)]
      false_positives = temp_sample[~temp_sample["tbprofiler_variant_substitution_nt"].isin(truth_mutations["tbprofiler_variant_substitution_nt"])]
      if not false_positives.empty:
        fp_sample = false_positives[false_positives["sample_id"].str.endswith(suffix)]

        fp_mutation = fp_sample["tbprofiler_variant_substitution_nt"].values[0]
        fp_gene_name = fp_sample["tbprofiler_gene_name"].values[0]
        fp_coverage_level = int(fp_sample["sample_id"].str.split("_").str[-1].values[0])
        fp_depth = int(fp_sample["depth"].values[0])
        fp_read_support = int(fp_sample["read_support"].values[0])
        fp_freq = float(fp_sample["frequency"].values[0])

        # Add lowest coverage FP subsample to the mutation_dictionary
        if fp_mutation not in mutation_dictionary:
          mutation_dictionary[fp_mutation] = [fp_mutation, fp_gene_name, fp_coverage_level, fp_depth, fp_read_support, fp_freq, "N/A", "N/A", "N/A"]
        elif fp_coverage_level < mutation_dictionary[fp_mutation][2]:
            mutation_dictionary[fp_mutation] = [fp_mutation, fp_gene_name, fp_coverage_level, fp_depth, fp_read_support, fp_freq, "N/A", "N/A", "N/A"]

    # The index of the last "False" value in this list will correspond to the highest coverage level where the mutation was NOT found.
    # Catches situations where the mutation disappears and then reappears - we want to count LOD at the coverage level right before it disappears for the first time.
    # If False is not found in the list, the mutation was not found in any subsamples (False Negative). Also calls false negatives if there are ZERO subsamples.
    try:
      last_false_index = len(true_positives) - true_positives[::-1].index(False) - 1

      # This is the subsample that has the lowest coverage level where the mutation was found (while considering disappearances).
      sample_choice = sorted_samples[last_false_index + 1]
    except (ValueError, IndexError) as e:
      mutation_dictionary[mutation] = [mutation, gene_name, "N/A", "N/A", "N/A", "N/A", truth_depth, truth_read_support, truth_freq]
      continue

    coverage_level = int(sample_choice.split("_")[-1])
    sample = one_sample[one_sample["sample_id"].str.endswith(sample_choice)]
    sample_depth = int(sample[sample["tbprofiler_variant_substitution_nt"] == mutation]["depth"].values[0])
    sample_freq = float(sample[sample["tbprofiler_variant_substitution_nt"] == mutation]["frequency"].values[0])
    sample_read_support = int(sample[sample["tbprofiler_variant_substitution_nt"] == mutation]["read_support"].values[0])
    mutation_dictionary[mutation] = [mutation, gene_name, coverage_level, sample_depth, sample_read_support, sample_freq, truth_depth, truth_read_support, truth_freq]
  return mutation_dictionary

def sort_df(original_df, prefix):
  original_df["coverage_level"] = original_df["sample_id"].str.split("_").str[-1].replace("full", 999).apply(lambda x: int(x))
  original_df["sample_name"] = original_df["sample_id"].str.split("_").str[0]
  #print(original_df["coverage_level"].unique())

  sorted_df = original_df.sort_values(by=["sample_name", "tbprofiler_variant_substitution_nt", "antimicrobial", "coverage_level"], ascending=[True, True, True, False], ignore_index=True)
  filename = prefix + ".sorted_laboratorian_reports.csv"
  sorted_df.to_csv(filename, index=False)

def main():
  options = get_input_args()
  # import the combined laboratorian reports
  original_df = pd.read_csv(options.laboratorian_reports, na_filter=False)

  # remove all WT mutations
  original_df = original_df[original_df["tbprofiler_variant_substitution_nt"] != "WT"]
  sort_df(original_df, options.output)

  suffixes = import_suffixes(options.suffix_file)

  suffixes.append(options.truth)
  # get list of sample IDs without any and all suffixes
  all_sample_ids = original_df["sample_id"]
  for suffix in suffixes:
    all_sample_ids = all_sample_ids.str.removesuffix(suffix)

  sample_ids = set(all_sample_ids.unique())
  #print("\n".join(sample_ids))
  # for each sample ID, extract the data and compare the suffixes
  lowest_coverage_levels = {}

  for sample in sample_ids:
    one_sample = original_df[original_df["sample_id"].str.startswith(sample)]
    # compare to the truth sample
    truth_sample = one_sample[one_sample["sample_id"].str.endswith(options.truth)]
    one_sample = one_sample[~one_sample["sample_id"].str.endswith(options.truth)]
    lowest_coverage_levels[sample] = compare_against_truth(one_sample, truth_sample)

  # sort by sample name first.
  lowest_coverage_levels = dict(sorted(lowest_coverage_levels.items()))

  # sort all the mutation_dictionaries by coverage level, placing FN and FP at the bottom.
  for sample in lowest_coverage_levels:
    lowest_coverage_levels[sample] = dict(sorted(lowest_coverage_levels[sample].items(), key=lambda item: item[1][2] if item[1][2] != "N/A" else float('inf')))

  # Loop through sorted dictionary and print out the results.
  with open("mutation_results.csv", "w") as f:
    header = ["Sample," "Mutation", "Gene," "Coverage Level," "Sample Depth," "Sample Read Support," "Sample Frequency," "Reference Depth," "Reference Read Support," "Reference Frequency"]
    print(",".join(header), file=f)

    for sample, results in lowest_coverage_levels.items():
      print(f"{sample}", file=f)
      print(f"Sample: {sample}")

      highest_coverage = max(results.items(), key=lambda item: item[1][2] if item[1][2] != "N/A" else float('inf'))
      highest_depth = max(results.items(), key=lambda item: item[1][3] if item[1][3] != "N/A" else float('inf'))
      highest_read_support = max(results.items(), key=lambda item: item[1][4] if item[1][4] != "N/A" else float('inf'))
      highest_freq = max(results.items(), key=lambda item: item[1][5] if item[1][5] != "N/A" else float('inf'))

      print(f"Highest coverage level where a mutation was found before FN result:\t {highest_coverage[1][2]} [{highest_coverage[1][0]}]")
      print(f"Highest sample depth where a mutation was found before FN result:\t {highest_depth[1][3]} [{highest_depth[1][0]}]")
      print(f"Highest read support where a mutation was found before FN result:\t {highest_read_support[1][4]} [{highest_read_support[1][0]}]")
      print(f"Highest allele frequency where a mutation was found before FN result:\t {highest_freq[1][5]} [{highest_freq[1][0]}]")

      for mutation, stats in results.items():
        extra_info = " ,"
        #False Negative
        if stats[2] == "N/A":
          extra_info = f"FN,"
          print(f"False Negative: {mutation} | Reference Depth: {stats[6]} | Reference Read Support: {stats[7]} | Reference Frequency: {stats[8]}")
        #False Positive
        elif stats[6] == "N/A":
          extra_info = f"FP,"
          print(f"False Positive: {mutation} | Coverage Level: {stats[2]} | Sample Depth: {stats[3]} | Sample Read Support: {stats[4]} | Sample Frequency: {stats[5]}")

        print(f"{extra_info}" + ",".join(str(_) for _ in stats), file=f)

      print("\n", file=f)
      print("\n")

if __name__ == "__main__":
  main()
