"""
TEST SCRIPT

Extract the mapping stats from STAR Log files. Useful for getting library sizes prior to differential expression analysis.
"""
import argparse
from glob import glob
from collections import defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('STAR_log_dir', help="The directory where the 'Log.final.out' files are saved")
    parser.add_argument('out_dir', help="Location to save the mapping stats.")
    args = parser.parse_args()

    star_dir = args.STAR_log_dir
    out_dir = args.out_dir

    star_files = glob(str(Path(star_dir, "*_Log.final.out")))
    mapping_stats = defaultdict(lambda : defaultdict(int))

    for star_file in star_files:
        sample_name = str(Path(star_file).stem).replace("_Log.final", "")
        with open(star_file, "r") as s:
            for line in s:
                l = line.strip()
                if l.startswith("Number of input reads"):
                    mapping_stats[sample_name]["input"] = l.split()[-1]
                elif l.startswith("Uniquely mapped reads number"):
                    mapping_stats[sample_name]["unique"] = l.split()[-1]
                elif l.startswith("Number of reads mapped to multiple loci"):
                    mapping_stats[sample_name]["multi"] = l.split()[-1]
                elif l.startswith("Number of reads unmapped: too short"):
                    mapping_stats[sample_name]["unmapped"] = l.split()[-1]

    if not Path(out_dir).exists():
        Path(out_dir).mkdir()

    outfile = Path(out_dir, "mapping_stats.tsv")
    with open(outfile, "w") as out:
        out.write(f"sample\tinput_reads\tuniquely_mapped\tmulti_mapped\tunmapped\n")
        for sample in mapping_stats.keys():
            out.write(f"{sample}\t{mapping_stats[sample]['input']}\t{mapping_stats[sample]['unique']}\t{mapping_stats[sample]['multi']}\t{mapping_stats[sample]['unmapped']}\n")


if __name__ == "__main__":
    main()