import subprocess
from collections import defaultdict
from pathlib import Path
import argparse


def compute_unique_coverage(unique_bam, repnames_bedfile):
    """Use bedtools to calculate the number of unique reads that overlap the repeat ranges"""
    print(f"Computing overlaps between {unique_bam} and {repnames_bedfile}...")
    parent_uniq_bam = Path(unique_bam).parent
    uniq_counts = Path(parent_uniq_bam, "_coverage.txt")
    command = f"bedtools coverage -counts -a {repnames_bedfile} -b {unique_bam} > {uniq_counts}"
    subprocess.run(command, shell=True)

    return Path(uniq_counts)


def count_unique(coverage_results, debug):
    """Create a dictionary of unique counts from the bedtools coverage output"""
    print("Computing unique counts...")
    counts = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
    with open(coverage_results, "r") as cov_res:
        for line in cov_res:
            l = line.strip().split()
            rep_element = l[3]
            class_ = l[4]
            family = l[5]
            count = int(l[6])
            
            counts[class_][family][rep_element] += count

    # remove the coverage results - unique counts output later in standard format
    if not debug:
        subprocess.run(f"rm {coverage_results}", shell=True)
    
    return counts


def count_multi(tmp_alignment_file, mapq):
    """Create a dictionary of counts for reads mapped to pseudogenome."""
    print("Counting reads overlapped to pseudogenome...")
    counts = defaultdict(float)

    with open(tmp_alignment_file, "r") as sam:
        for line in sam:
            if line.startswith("@"):
                continue  # skip headers
            else:
                l = line.strip().split("\t")
                rep_element = l[2]
                flag = int(l[1])
                map_q = int(l[4])
                if map_q >= mapq:
                    counts[rep_element] += 1

    return counts


def count_fractional(tmp_alignment_file, mapq):
    """Return the fractional counts from the pseudo aligned reads by dividing the count of the elements
    by the number of subfamilies the read aligns to."""
    print("Computing fractional counts...")
    counts = defaultdict(lambda : defaultdict(float))
    frac_counts = defaultdict(float)

    with open(tmp_alignment_file, "r") as sam:
        for line in sam:
            if line.startswith("@"):
                continue  # skip headers
            else:
                l = line.strip().split("\t")
                read_id = l[0]
                flag = int(l[1])
                rep_element = l[2]
                map_q = int(l[4])
                if map_q >= mapq: 
                    counts[read_id][rep_element] = 1  # every element gets numerator 1 (even if observed multiple times)

    for read_id in counts.keys():
        N_subfamilies = len(counts[read_id].keys())
        for sub_family in counts[read_id].keys():
            sub_family_count = counts[read_id][sub_family]  # count should be 1 as defined above
            frac_counts[sub_family] += sub_family_count / N_subfamilies  # frac count == 1 / N_subfamilies read aligns to

    return frac_counts


def create_element_mapping(repnames_bedfile):
    """Create a mapping of the element names to their classes and families"""
    elem_key = defaultdict(lambda : defaultdict(str))

    with open(repnames_bedfile, "r") as bed:
        for line in bed:
            l = line.strip().split("\t")
            name = l[3]
            class_ = l[4]
            family = l[5]
            elem_key[name]["class"] = class_
            elem_key[name]["family"] = family
    
    return elem_key  


def map_counts(multi_counts, element_mapper):
    """Create a dictionary of multi-mapped reads using the class > family > element hierarchy.
    Since the output from the count_multi function returns element names and their counts, 
    these element names must be re-associated with the classes and families they are part of. This 
    association is performed by the element_mapper.
    """
    hierarchy = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
    for elem in multi_counts.keys():
        class_ = element_mapper[elem]['class']
        family = element_mapper[elem]['family']
        hierarchy[class_][family][elem] += multi_counts[elem]

    return hierarchy


def align_to_pseudogenome(pseudogenome_STAR_idx, multi_fastq, element_mapper, threads, debug, mapq):
    """Align the multi mapping reads to the pseudo genome and return counts of reads uniquely mapping to pseudogenome"""
    multi_fastq_parent = Path(multi_fastq).parent
    out_dir = Path(multi_fastq_parent, "__tmp_pseudo_aligned")
    out_prefix = Path(out_dir, "__tmp_")

    if not Path(out_dir).exists():
        Path(out_dir).mkdir()
    
    print(f"Aligning {multi_fastq} to pseudogenome using STAR...")
    align_cmd = f"STAR --runThreadN {threads} --genomeDir {pseudogenome_STAR_idx} --readFilesIn {multi_fastq} --outFileNamePrefix {out_prefix}"
    subprocess.run(align_cmd, shell=True)

    tmp_alignment_file = Path(out_dir, "__tmp_Aligned.out.sam")

    # first, count all of the multi-mapped reads
    counts = count_multi(tmp_alignment_file, mapq)
    mapped_multi_counts = map_counts(counts, element_mapper)

    # then use the same alignment file to find the fractional counts
    frac_counts = count_fractional(tmp_alignment_file, mapq)
    mapped_frac_counts = map_counts(frac_counts, element_mapper)

    # clean the temporary alignment directory
    if not debug:
        subprocess.run(f"rm -r {out_dir}", shell=True)

    return mapped_multi_counts, mapped_frac_counts


def combine_counts(counts1, counts2):
    """Generic function for combining two dictionaries of count data"""
    total_counts = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))

    for class_ in counts1.keys():
        for family in counts1[class_].keys():
            for elem in counts1[class_][family].keys():
                total_counts[class_][family][elem] += counts1[class_][family][elem]
        
    for class_ in counts2.keys():
        for family in counts2[class_].keys():
            for elem in counts2[class_][family].keys():
                total_counts[class_][family][elem] += counts2[class_][family][elem]
            
    return total_counts


def write_output(count_data, outfile):
    """Write count data to the outfile. Output truncates floats to ints."""
    print(f"Writing results to {outfile}...")
    with open(outfile, "w") as f:
        for class_ in count_data:
            for family in count_data[class_]:
                for elem in count_data[class_][family]:
                    count = int(count_data[class_][family][elem])
                    f.write(f"{class_}\t{family}\t{elem}\t{count}\n")


def summarize(count_data):
    """Collapse the counts to the class and family level."""
    class_data = defaultdict(float)
    family_data = defaultdict(lambda : defaultdict(float))
    
    # class level counts
    for class_ in count_data:
        for family in count_data[class_]:
            for elem in count_data[class_][family]:
                count = int(count_data[class_][family][elem])
                class_data[class_] += count
    
    # family level counts
    for class_ in count_data:
        for family in count_data[class_]:
            for elem in count_data[class_][family]:
                count = int(count_data[class_][family][elem])
                family_data[class_][family] += count
    
    return class_data, family_data


def write_class_summary(class_counts, outfile):
    """Write class-level counts to outfile. Output truncates floats to ints"""
    with open(outfile, "w") as f:
        for class_ in class_counts.keys():
            count = int(class_counts[class_])
            f.write(f"{class_}\t{count}\n")


def write_family_summary(family_counts, outfile):
    """Write family-level counts to outfile. Output truncates floats to ints."""
    with open(outfile, "w") as f:
        for class_ in family_counts.keys():
            for family in family_counts[class_]:
                count = int(family_counts[class_][family])
                f.write(f"{class_}\t{family}\t{count}\n")


def main():
    parser = argparse.ArgumentParser(description='Quantifying repeat counts by mapping multimapped reads to pseudogenomes and counting overlaps with repeat elements for uniquely mapped reads.',
    usage = 'python fastRE_count.py sampleName path/to/fastRE_Setup path/to/sampleName_unique.bam --threads 8 --summarize')
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('sampleName',help='The name of the sample to be processed. Typically the prefix of the unique.bam file.')
    parser.add_argument('setupFolder',help='path to the fastRE_Setup folder')
    parser.add_argument('uniqueAlignmentFile', help='path/to/sampleName_unique.bam')
    parser.add_argument('--mapq', default=0, metavar=0, type=int, help="""The MAPQ threshold for multimapped reads, Should be one of [255, 3, 2, 1, 0]. 
    STAR uses 255 to designate uniquely mapped reads. 3=Read maps to 2 locations, 2=Read maps to 3 locations, 1=Read maps to 4-9 locations, 0=Read maps 
    to 10 or more locations. The default is set to 0 to allow for all multimapped that align to the pseudogenome to be counted. This setting diverges 
    from the original RepEnrich2 program since the aligners (STAR vs. bowtie2) use different MAPQ schemes.""")
    parser.add_argument('--pairedEnd', dest='pairedEnd', action='store_true', help='Designate this option for paired-end sequencing.')
    parser.add_argument('--threads', default=1, type=int, metavar=1, help='Number of threads to use for alignment with STAR.')
    parser.add_argument('--summarize', dest='summarize', action='store_true', help='In addition to the repeat name level output, produce collapsed counts at the class and family levels for each of the count types.')
    parser.add_argument('--debug', dest='debug', action='store_true', help='Select this option to prevent the removal of temporary files; useful for debugging')
    parser.set_defaults(pairedEnd=False, summarize=False, debug=False)
    args = parser.parse_args()

    sample_name = args.sampleName
    setup_folder = args.setupFolder
    alignment_file = args.uniqueAlignmentFile
    paired_end = args.pairedEnd
    threads = args.threads
    debug = args.debug
    mapq = args.mapq

    # main count routine ------------------------------------------------------
    # setup outfile paths - saved in same directory as the alignment file
    out_parent = Path(alignment_file).parent
    repnames_bedfile = Path(setup_folder, "repnames.bed")
    pseudogenome_STAR_idx = Path(setup_folder, "STAR_pseudogenome_idx")
    se_fastq = Path(out_parent, sample_name + "_multimap.fastq")   # The multi fastq files are expected to be in the same directory
    fastq1 = Path(out_parent, sample_name + "_multimap_R1.fastq")  #   as the unique.bam file. This spec may change in the future
    fastq2 = Path(out_parent, sample_name + "_multimap_R2.fastq")
    total_counts_outfile = Path(out_parent, sample_name + "_total_counts.tsv")
    uniq_counts_outfile = Path(out_parent, sample_name + "_unique_counts.tsv")
    fractional_counts_outfile = Path(out_parent, sample_name + "_fractional_counts.tsv")
    class_total_counts_outfile = Path(out_parent, sample_name + "_class_total_counts.tsv")
    class_uniq_counts_outfile = Path(out_parent, sample_name + "_class_unique_counts.tsv")
    class_fractional_counts_outfile = Path(out_parent, sample_name + "_class_fractional_counts.tsv")
    family_total_counts_outfile = Path(out_parent, sample_name + "_family_total_counts.tsv")
    family_uniq_counts_outfile = Path(out_parent, sample_name + "_family_unique_counts.tsv")
    family_fractional_counts_outfile = Path(out_parent, sample_name + "_family_fractional_counts.tsv")

    # run the routine
    print(f"Begin processing {sample_name}...")
    cov_res = compute_unique_coverage(alignment_file, repnames_bedfile)
    uniq_counts = count_unique(cov_res, debug)
    elem_mapper = create_element_mapping(repnames_bedfile)

    if paired_end:
        multi_count1, frac_count1 = align_to_pseudogenome(pseudogenome_STAR_idx, fastq1, elem_mapper, threads, debug, mapq)
        multi_count2, frac_count2 = align_to_pseudogenome(pseudogenome_STAR_idx, fastq2, elem_mapper, threads, debug, mapq)
        multi_counts = combine_counts(multi_count1, multi_count2)
        frac_counts = combine_counts(frac_count1, frac_count2)
    else:
        multi_counts, frac_counts = align_to_pseudogenome(pseudogenome_STAR_idx, se_fastq, elem_mapper, threads, debug, mapq)

    total_counts = combine_counts(uniq_counts, multi_counts)
    fractional_counts = combine_counts(uniq_counts, frac_counts)
    write_output(fractional_counts, fractional_counts_outfile)
    write_output(total_counts, total_counts_outfile)
    write_output(uniq_counts, uniq_counts_outfile)

    if summarize:
        print("Writing class-level and family-level summarized counts...")
        class_total, family_total = summarize(total_counts)
        class_uniq, family_uniq = summarize(uniq_counts)
        class_fractional, family_fractional = summarize(fractional_counts)
        write_class_summary(class_total, class_total_counts_outfile)
        write_class_summary(class_uniq, class_uniq_counts_outfile)
        write_class_summary(class_fractional, class_fractional_counts_outfile)
        write_family_summary(family_total, family_total_counts_outfile)
        write_family_summary(family_uniq, family_uniq_counts_outfile)
        write_family_summary(family_fractional, family_fractional_counts_outfile)

    print("Finished.")
    print("#" * 80)


if __name__ == '__main__':
    main()