"""
TEST SCRIPT

Count reads aligned to consensus without the bedfile annotation
"""
from collections import defaultdict
from pathlib import Path
import argparse

import pysam


def combine_counts(count_data1, count_data2):
    """Combine unique and frac counts for final fractional counting count_data"""
    total = defaultdict(float)
    for subfamily in count_data1:
        total[subfamily] += count_data1[subfamily]
    
    for subfamily in count_data2:
        total[subfamily] += count_data2[subfamily]

    return total


def count(bamfile):
    """Return Unique, fractional and total counts from parsing the bamfile"""
    print("Counting reads overlapped to pseudogenome...")
    uniq_counts = defaultdict(float)
    total_counts = defaultdict(float)
    frac_numerator = defaultdict(lambda: defaultdict(float))
    frac_counts = defaultdict(float)

    # calculate uniq counts and tally multi counts
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for segment in bam:
            mapq = segment.mapping_quality
            rep_element = segment.reference_name
            read_id = segment.query_name

            # tally total counts for all rep elemnets in the alignment
            total_counts[rep_element] += 1

            if mapq == 255:  # calculate unique counts
                uniq_counts[rep_element] += 1
            else:            # tally the sub families that a read multimaps to
                frac_numerator[read_id][rep_element] = 1  # every element gets numerator 1 (even if observed multiple times)

    # calculating fractional counts
    for read_id in frac_numerator.keys():
        N_subfamilies = len(frac_numerator[read_id].keys())
        for sub_family in frac_numerator[read_id].keys():
            frac_counts[sub_family] += (1 / N_subfamilies)  # frac count == 1 / N_subfamilies read aligns to

    return uniq_counts, frac_counts, total_counts


def write_output(count_data, outfile):
    """Write count data to the outfile. Output truncates floats to ints."""
    print(f"Writing results to {outfile}...")
    with open(outfile, "w") as f:
        for subfamily in count_data:
            count = int(count_data[subfamily])
            f.write(f"{subfamily}\t{count}\n")


def main():
    parser = argparse.ArgumentParser(prog="count_from_pseudoalignment.py",
                                    usage="count_from_pseudoalignment.py file_prefix in.sam repnames.bed")
    parser.add_argument('file_prefix', help="Sample name to prefix the results files with")
    parser.add_argument('in_bam', help="SAM file of reads aligned to the pseudogenome")
    parser.add_argument('--out_dir', help="Location to save results. Defaults to parent directory of the bam file.")
    args = parser.parse_args()

    file_prefix = args.file_prefix
    bamfile = args.in_bam
    out_dir = args.out_dir

    # set up outfile locations ------------------------------------------------
    if out_dir:
        out_directory = Path(out_dir)
        if not out_directory.exists():
            out_directory.mkdir()
    else:
        out_directory = Path(bamfile).parent

    uniq_out = Path(out_directory, f"{file_prefix}_unique_counts.tsv")
    frac_out = Path(out_directory, f"{file_prefix}_fractional_counts.tsv")
    total_out = Path(out_directory, f"{file_prefix}_total_counts.tsv")

    # perform counting and write output ---------------------------------------
    uniq, frac, tot = count(bamfile)
    fractional_counts = combine_counts(uniq, frac)
    
    write_output(uniq, uniq_out)
    write_output(fractional_counts, frac_out)
    write_output(tot, total_out)


if __name__ == '__main__':
    main()
