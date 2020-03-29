import argparse
from shlex import split
import subprocess
from pathlib import Path
from os import devnull, getcwd

import pysam
from pybedtools import BedTool


def get_unique(sample_name, alignment_file, threads, mapq, bwt_mode, paired_end):
    """Subset unique reads from alignment file"""
    algn_parent = Path(alignment_file).parent
    uniq = Path(algn_parent, sample_name + "_unique.bam")

    if bwt_mode:
        if paired_end:
            pysam.view("-f", "3", "-q", f"{mapq}", "-b", "-o", f"{uniq}", f"{alignment_file}", catch_stdout=False)
        else:
            pysam.view("-F", "4", "-q", f"{mapq}", "-b", "-o", f"{uniq}", f"{alignment_file}", catch_stdout=False)
    else:
        pysam.view("--threads", f"{threads}", "-q", f"{mapq}", "-b", "-o", f"{uniq}", f"{alignment_file}", catch_stdout=False)

    print(f"Subsetting unique reads from {alignment_file}...")
    if bwt_mode:
        print(f"\tBowtie2 mode selected. Using {mapq} as the MAPQ quality score to filter on.")


def get_multimapped(alignment_file, threads, mapq, bwt_mode, paired_end):
    """Subset multi-mapped reads from the alignment file"""
    algn_parent = Path(alignment_file).parent
    multi_bam = Path(algn_parent, "_multi.bam")
    proper_pair = Path(algn_parent, "_proper_pair.bam")

    print(f"Subsetting Multi-mapped reads...")
    if bwt_mode:
        if paired_end:
            pysam.view("--threads", f"{threads}", "-f", "3", "-b", "-o", f"{proper_pair}", f"{alignment_file}", catch_stdout=False)
            pysam.view("--threads", f"{threads}", "-q", f"{mapq}", "-b", "-U", f"{multi_bam}", f"{proper_pair}", catch_stdout=True)
        else:
            pysam.view("--threads", f"{threads}", "-F", "4", "-q", f"{mapq}", "-b", "-U", f"{multi_bam}", f"{alignment_file}", catch_stdout=True)
    else:
        pysam.view("--threads", f"{threads}", "-q", f"{mapq}", "-b", "-U", f"{multi_bam}", f"{alignment_file}", catch_stdout=True)

    return multi_bam


def multi_to_fastq(sample_name, multi_bam, paired_end, alignment_file, threads):
    algn_parent = Path(alignment_file).parent
    srt_multi_bam = Path(algn_parent, "_sorted.multi.bam")
    out_SE = Path(algn_parent, sample_name + "_multimap.fastq")
    out_R1 = Path(algn_parent, sample_name + "_multimap_R1.fastq")
    out_R2 = Path(algn_parent, sample_name + "_multimap_R2.fastq")
    
    if not paired_end:
        print(f"Writing multi-mapped reads to {out_SE}...")
        bamfile = BedTool(multi_bam)
        BedTool.bam_to_fastq(bamfile, fq=out_SE)
    else:
        print("Sorting multi-mapped reads...")
        pysam.sort("--threads", f"{threads}", "-n", "-o", f"{srt_multi_bam}", f"{multi_bam}")
        sorted_bamfile = BedTool(srt_multi_bam)
        print(f"Writing multi-mapped reads to {out_R1} {out_R2}...")
        BedTool.bam_to_fastq(sorted_bamfile, fq=out_R1, fq2=out_R2)


def cleanup(alignment_file, debug):
    """remove intermediate files"""
    algn_parent = Path(alignment_file).parent
    srt_multi_bam = Path(algn_parent, "_sorted.multi.bam")
    multi_bam = Path(algn_parent, "_multi.bam")
    pp_bam = Path(algn_parent, "_proper_pair.bam")

    if not debug:
        print("Performing cleanup of large intermediate files...")
        subprocess.run(f"rm -f {multi_bam} {srt_multi_bam} {pp_bam}", shell=True)


def main():
    parser = argparse.ArgumentParser(description="Prepartion for downstream analysis with fastRE_count.py. Subsets unique and multi-mapped reads and creates fastq files for multi-mapped reads.", 
    usage = "python fastRE_subset.py path/to/sampleName_Aligned.out.bam sampleName --threads 8")
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('sampleName', help='The name of your sample (used to generate the name of the output files)')
    parser.add_argument('alignmentFile',help='bam file of aligned reads. Example : path/to/sampleName_Aligned.out.bam')
    parser.add_argument('--pairedEnd', dest='pairedEnd', action='store_true', help='Designate this option for paired-end data.')
    parser.add_argument('--threads', default=0, type=int, metavar=0, help='Additional number of threads to use in samtools calls.')
    parser.add_argument('--debug', dest='debug', action='store_true', help='Select this option to prevent the removal of temporary files; useful for debugging')
    parser.add_argument('--bowtieMode', action='store_true', help="Set this flag if you would like to use bowtie2 instead of STAR for all downstream analyses.")
    parser.add_argument('--MAPQ', default=255, type=int, metavar=255, help="Set the MAPQ score for uniquely mapping reads. In STAR this value is 255. If bowtieMode is set then this value MUST be changed (30 is the default in RepEnrich2).")
    parser.set_defaults(pairedEnd=False, debug=False, bowtieMode=False)
    args = parser.parse_args()

    alignment_file = args.alignmentFile
    sample_name = args.sampleName
    paired_end = args.pairedEnd
    threads = args.threads
    debug = args.debug
    bwt_mode = args.bowtieMode
    mapq = args.MAPQ
    
    # main subset routine -----------------------------------------------------
    get_unique(sample_name, alignment_file, threads, mapq, bwt_mode, paired_end)
    multi_path = get_multimapped(alignment_file, threads, mapq, bwt_mode, paired_end)
    multi_to_fastq(sample_name, multi_path, paired_end, alignment_file, threads)
    cleanup(alignment_file, debug)
    print(f"Finished processing {sample_name}.")
    print("#" * 80)


if __name__ == '__main__':
    main()