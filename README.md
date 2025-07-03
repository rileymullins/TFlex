# multiplexed_srt_tf_mapping
Identify peaks from multiplexed self-reporting transposon data for identifying transcription factor binding sites.

# Annotation file format
This file should be a .csv or .tsv format. The columns must be in this order: library_name, sample_barcode, sample_name, group_name. Here is an example tab-separated structure:

library_name	sample_barcode	sample_name	group_name
Library_A	AAGGCAGACG	Rep1_HyPBase  HyPBase
Library_A	TACCGCTGAC	Rep1_TF_A	TF_A
Library_A	AAGATTAGAC	Rep1_TF_B	TF_B
Library_A	AGTATGACCG	Rep1_TF_C	TF_C
Library_A	TACTAGGTTG	Rep1_TF_D	TF_D
Library_A	CTCTGGATTA	Rep1_TF_E	TF_E
Library_A	AGTTAACGAT	Rep1_TF_F	TF_F
Library_A	GAGACACGTC	Rep1_TF_G	TF_G
Library_A	GTAATGCTTC	Rep1_TF_H	TF_H
Library_A	TTAACGATCG 	Rep1_TF_I	TF_I
Library_B	AAGGCAGACG	Rep2_HyPBase	HyPBase
Library_B	TACCGCTGAC	Rep2_TF_A	TF_A
Library_B	AAGATTAGAC	Rep2_TF_B	TF_B
Library_B	AGTATGACCG	Rep2_TF_C	TF_C
Library_B	TACTAGGTTG	Rep2_TF_D	TF_D
Library_B	CTCTGGATTA	Rep2_TF_E	TF_E
Library_B	AGTTAACGAT	Rep2_TF_F	TF_F
Library_B	GAGACACGTC	Rep2_TF_G	TF_G
Library_B	GTAATGCTTC	Rep2_TF_H	TF_H
Library_B	TTAACGATCG 	Rep2_TF_I	TF_I


# Variable names
consensus_peak_id = chrom:consensus_peak_start-consensus_peak_end  
chrom = chromosome  
consensus_peak_start = start coordinate of final merged, pan-dataset consensus peak  
consensus_peak_end = end coordinate of final merged, pan-dataset consensus peak  
consensus_peak_width = width (end – start) of final merged, pan-dataset consensus peak  

num_samples_in_consensus_peak = number of unique sample_name values in this consensus peak  
num_groups_in_consensus_peak = number of unique group_name values in this consensus peak  

sample_name = sample identifier (from annotation file). All sample_names MUST BE UNIQUE. At the most basic level, this identifier is the peak calling unit on which peaks will be called independently. This value could be a replicate of the same TF (replicate-1_TF-A), a unique time point of the same TF (timepoint1_replicate1_TF-A), or an otherwise unique experimental condition (drugX_replicate1_TF-A).

group_name = group identifier (from annotation file). At the most basic level, this identifier is used to aggregate stats across samples belong to the same broader group. Matching the examples above, the group_name for a replicate of the same TF may be "TF_A", group_name for a time point of the same TF may be "timepoint1_TF-A", and group_nae for an otherwise unique experimental condition may be "drugX_TF-A")

fragment_peak_start_for_sample = fragment-based peak start coordinate for this sample_name  
fragment_peak_end_for_sample = fragment-based peak end coordinate for this sample_name  
fragment_peak_width_for_sample = fragment-based peak width (end – start) for this sample_name  

srt_bc_peak_start_for_sample = SRT-barcode-based peak start coordinate for this sample_name  
srt_bc_peak_end_for_sample = SRT-barcode-based peak end coordinate for this sample_name  
srt_bc_peak_width_for_sample = SRT-barcode-based peak width (end – start) for this sample_name  

sample_total_reads_in_consensus_peak = total read count in consensus peak for this sample_name  
sample_total_fragments_in_consensus_peak = total fragment count (merged qbed rows after SRT-barcode correction and deduplication) in consensus peak for this sample_name  
sample_total_insertions_in_consensus_peak = total unique insertions (unique SRT barcodes) in consensus peak for this sample_name  

group_total_reads_in_consensus_peak = sum of sample_total_reads_in_consensus_peak across all sample_name values in this group_name  
group_total_fragments_in_consensus_peak = sum of sample_total_fragments_in_consensus_peak across all sample_name values in this group_name  
group_total_insertions_in_consensus_peak = sum of sample_total_insertions_in_consensus_peak across all sample_name values in this group_name  


# Input QBED format
Column 1 = chromosome
Column 2 = Start
Colu






