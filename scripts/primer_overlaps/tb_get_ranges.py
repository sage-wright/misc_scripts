import pandas as pd

primer_file = pd.read_csv("/home/sage_wright/github/misc_scripts/scripts/primer_overlaps/tb_primers-nymain.bed", sep="\t", header=None)
primers = pd.DataFrame(primer_file)

# extract regions from primer file
regions = primers.iloc[:,1:4]
regions.columns = ['Forward Start', 'Reverse End', 'Gene']
regions.set_index('Gene', inplace=True)

regions_dict = regions.T.to_dict('list')

# dictionary : {target: [forward (F) start, reverse (R) end]}
ranges_dict = {}

for target, region in regions_dict.items():
  # match forward targets (end with F) with reverse targets (end with R)
  # in order to make dictionary of ranges
  if target.endswith('F'):
    # get target name by removing suffix (everything after _ or -)
    name = target.rsplit('-', 1)[0]

    # find matching name with R suffix
    try:
      forward_start = region[0]
      reverse_end = regions_dict[name + '-R'][1]

      ranges_dict[name] = [str(forward_start), str(reverse_end)]    
    except:
      ranges_dict[name] = [region[0], None]

for target, region in ranges_dict.items():
  print(target + '\t' + '\t'.join(region))
