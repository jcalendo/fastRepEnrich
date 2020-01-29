import argparse
from shlex import split
import subprocess
from pathlib import Path
from os import devnull, getcwd


def get_unique(sample_name, alignment_file, threads):
    """Subset unique reads from alignment file"""
    algn_parent = Path(alignment_file).parent
    uniq = Path(algn_parent, sample_name + "_unique.bam")
    uniq_cmd = f"samtools view --threads {threads} -q 255 -b -o {uniq} {alignment_file}"

    print(f"Subsetting unique reads from {alignment_file}...")
    subprocess.run(uniq_cmd, shell=True)


def get_multimapped(alignment_file, threads):
    """Subset multi-mapped reads from the alignment file"""
    algn_parent = Path(alignment_file).parent
    multi_bam = Path(algn_parent, "_multi.bam")
    multi_cmd = f"samtools view --threads {threads} -q 255 -b -U {multi_bam} {alignment_file} > {devnull}"
    
    print(f"Subsetting Multi-mapped reads...")
    subprocess.run(multi_cmd, shell=True)

    return multi_bam


def multi_to_fastq(sample_name, multi_bam, paired_end, alignment_file, threads):
    algn_parent = Path(alignment_file).parent
    srt_multi_bam = Path(algn_parent, "_sorted.multi.bam")
    out_SE = Path(algn_parent, sample_name + "_multimap.fastq")
    out_R1 = Path(algn_parent, sample_name + "_multimap_R1.fastq")
    out_R2 = Path(algn_parent, sample_name + "_multimap_R2.fastq")
    if not paired_end:
        print(f"Writing multi-mapped reads to {out_SE}...")
        fastq_cmd = f"samtools fastq --threads {threads} -F 0x900 {multi_bam} > {out_SE}"
    else:
        print("Sorting multi-mapped reads...")
        srt_cmd = f"samtools sort --threads {threads} -n -o {srt_multi_bam} {multi_bam}"
        subprocess.run(srt_cmd, shell=True)
        print(f"Writing multi-mapped reads to {out_R1} {out_R2}...")
        fastq_cmd = f"samtools fastq --threads {threads} -1 {out_R1} -2 {out_R2} -0 {devnull} -s {devnull} -n -F 0x900 {srt_multi_bam}"
    
    subprocess.run(fastq_cmd, shell=True)


def cleanup(alignment_file, debug):
    """remove intermediate files"""
    algn_parent = Path(alignment_file).parent
    srt_multi_bam = Path(algn_parent, "_sorted.multi.bam")
    multi_bam = Path(algn_parent, "_multi.bam")

    if not debug:
        print("Performing cleanup of large intermediate files...")
        subprocess.run(f"rm -f {multi_bam} {srt_multi_bam}", shell=True)


def main():
    parser = argparse.ArgumentParser(description="Prepartion for downstream analysis with fastRE_count.py. Subsets unique and multi-mapped reads and creates fastq files for multi-mapped reads.", 
    usage = "python fastRE_subset.py path/to/sampleName_Aligned.out.bam sampleName --threads 8")
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('alignmentFile',help='bam file of aligned reads. Example : path/to/sampleName_Aligned.out.bam')
    parser.add_argument('sampleName', help='The name of your sample (used to generate the name of the output files)')
    parser.add_argument('--pairedEnd', dest='pairedEnd', action='store_true', help='Designate this option for paired-end data.')
    parser.add_argument('--threads', default=0, type=int, metavar=0, help='Additional number of threads to use in samtools calls.')
    parser.add_argument('--debug', dest='debug', action='store_true', help='Select this option to prevent the removal of temporary files; useful for debugging')
    parser.set_defaults(pairedEnd=False, debug=False)
    args = parser.parse_args()

    alignment_file = args.alignmentFile
    sample_name = args.sampleName
    paired_end = args.pairedEnd
    threads = args.threads
    debug = args.debug

    # main subset routine -----------------------------------------------------
    get_unique(sample_name, alignment_file, threads)
    multi_path = get_multimapped(alignment_file, threads)
    multi_to_fastq(sample_name, multi_path, paired_end, alignment_file, threads)
    cleanup(alignment_file, debug)
    print(f"Finished processing {sample_name}.")
    print("#" * 80)


if __name__ == '__main__':
    main()