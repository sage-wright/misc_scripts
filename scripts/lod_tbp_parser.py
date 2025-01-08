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
  truth_mutations = truth_sample[["tbprofiler_variant_substitution_nt", "read_support"]]
  truth_mutations = truth_mutations.drop_duplicates(subset=['tbprofiler_variant_substitution_nt'], ignore_index=True)
    
  # this dictionary has the mutation as the key and the HIGHEST coverage level where the mutation was NOT found and the mutation's read support at that coverage level
  mutation_dictionary = {}
  
  for mutation in truth_mutations["tbprofiler_variant_substitution_nt"]:
    read_support = truth_mutations[truth_mutations["tbprofiler_variant_substitution_nt"] == mutation]["read_support"].values[0]

    for suffix in one_sample["sample_id"].unique():
      if mutation not in one_sample[one_sample["sample_id"].str.endswith(suffix)]["tbprofiler_variant_substitution_nt"].unique():
        # determine the highest coverage level where the mutation is NOT found (the last integer in the suffix [_#])        
        coverage_level = int(suffix.split("_")[-1])

        if mutation in mutation_dictionary:
          if coverage_level > mutation_dictionary[mutation][0]:
            mutation_dictionary[mutation] = [coverage_level, read_support]
        else:
          mutation_dictionary[mutation] = [coverage_level, read_support]
  
  
  print(truth_sample["sample_id"].unique()[0], ":", mutation_dictionary)  
  try:
    print("Highest coverage level where a mutation was NOT found:", max(mutation_dictionary.values()))
    return max(mutation_dictionary.values())

  except:
    if (len(truth_mutations["tbprofiler_variant_substitution_nt"]) == len(mutation_dictionary)):
      print("the number of truth mutations ({}) matched the number of mutations in the mutation_dictionary ({})".format(len(truth_mutations["tbprofiler_variant_substitution_nt"]), len(mutation_dictionary)))
      return 0
    else:
      print("CAUTION!!!!!!! the number of truth mutations ({}) did not match the number of mutations in the mutation_dictionary ({})".format(len(truth_mutations["tbprofiler_variant_substitution_nt"]), len(mutation_dictionary)))
      print("these mutations were not found in the subsampled data:")
      print(truth_mutations)
      return 999
    
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
  
  highest_coverage_levels = {}
  
  print("Hello! You're running `lod_tbp_parser.py`\n")
  print("the following lists take the following format: ")
  print("[reference sample name] {'mutation': [coverage level where mutation was not found, 'read support for that mutation in the reference sample'], ...}\n")
  
  for sample in sample_ids:
    one_sample = original_df[original_df["sample_id"].str.startswith(sample)]

    # compare to the truth sample
    truth_sample = one_sample[one_sample["sample_id"].str.endswith(options.truth)]
    one_sample = one_sample[~one_sample["sample_id"].str.endswith(options.truth)]
    highest_coverage_levels[sample] = compare_against_truth(one_sample, truth_sample)
    print()
      
  print("\nthe following dictionary shows the coverage level where drop-out occurs for each sample\n")
  print("it has the following format:")
  print("sample: [coverage level where drop-out occurs, read support at that coverage level]")
  print("OR")
  print("sample: coverage level where drop-out occurs (999 indicates certain mutations were only found in the reference sample)\n")
  print(highest_coverage_levels)
  print("coverage level where drop-out occurs: " + str(max(result[0] if type(result) is list else result for result in highest_coverage_levels.values())))

if __name__ == "__main__":
  main()
