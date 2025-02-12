# format: gene: min, max
overlaps = {}

# where a < b, and c < d
def check_overlap(a, b, c, d):
  return max(a, c) <= min(b, d)

with open("/home/sage_wright/github/misc_scripts/scripts/primer_overlaps/tb-ranges-nymain.tsv") as ranges:
  for line in ranges:
    line = line.strip().split('\t')
    gene = line[0]
    forward_start = int(line[1])
    try:
      reverse_end = int(line[2])
    except:
      overlaps[gene + "_2"] = [forward_start, None]
      continue
    
    if gene not in overlaps.keys():
      min_start = min(forward_start, reverse_end)
      max_end = max(forward_start, reverse_end)
      overlaps[gene] = [min_start, max_end]
    else:
      region_min_start = min(forward_start, reverse_end)
      region_max_end = max(forward_start, reverse_end)
  
      if check_overlap(region_min_start, region_max_end, overlaps[gene][0], overlaps[gene][1]):
        # update min and max values
        overlaps[gene][0] = min(region_min_start, overlaps[gene][0])
        overlaps[gene][1] = max(region_max_end, overlaps[gene][1])
      
      else:
        # make new entry for gene
        overlaps[gene + "_2"] = [region_min_start, region_max_end]
      
for gene, region in overlaps.items():
  gene = gene.strip()
  start = str(region[0])
  end = str(region[1])
  print(gene + '\t' + start + '\t' + end)