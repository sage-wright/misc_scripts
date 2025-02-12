# Coverage Locations

This directory contains files relevant to a coverage analysis I performed to test the efficacy of various primers for M.tb.

Two scripts are kept here:

1. `get_full_primer_overlap_range.py` is used on the modified output of the scripts contained in ../scripts/primer_overlaps/ to extract the full range that the primer covers.

    - the input tab-delimited file takes the following format (no header):
      - `primer_name	start_genomic_position	end_genomic_position	target_gene`
    - this will generate a tab-delimited file with the following format (no header):
      - `primer_name	genomic_position	target_gene`

      for this example, the primer name and target gene are repeated for every genomic_position where that primer maps (i.e., every position in the range(start_genomic_position, end_genomic_position))

2. `coverage_assessment.py` is used to calculate regions of > 20x depth for the target regions; in this case, the primers's targeted region. It will also generate figures for each gene target(see the `example_high-coverage_ranges.png` file that is included in this directory as an example).

    - In order to run this script, the output from `get_full_primer_overlap_range.py` must be slightly modified to include the tab-delimited header:
      - `BAM File	Position	Gene`

      In addition, any forward or reverse tags should be removed (e.g., `primerA-F` and `primerA-R` should _both_ turn into `primerA`) so that the primers appear on the same line.

    - An additional input is the file produced in the `scripts/primer_overlaps` directory after running the two included scripts. Please see the README in that directory for more details.

  The output of this task is a series of printed statements and a figure for each gene target that shows the regions of > 20x depth. The figure is saved as a .png file in the same directory where the script is run.

  It is recommended to redirect the stdout to a file.

Finally, there is the H37Rv.gb file for reference.
  