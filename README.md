# multiplexed_srt_tf_mapping
Identify peaks from multiplexed self-reporting transposon data for identifying transcription factor binding sites.

# Variable names
consensus_peak_id = chrom:consensus_peak_start-consensus_peak_end  
chrom = chromosome  
consensus_peak_start = start coordinate of final merged, pan-dataset consensus peak  
consensus_peak_end = end coordinate of final merged, pan-dataset consensus peak  
consensus_peak_width = width (end – start) of final merged, pan-dataset consensus peak  

num_samples_in_consensus_peak = number of unique sample_name values in this consensus peak  
num_groups_in_consensus_peak = number of unique group_name values in this consensus peak  

sample_name = sample identifier (from annotation file)  
group_name = group identifier (from annotation file)  

frag_peak_start_for_sample = fragment-based peak start coordinate for this sample_name  
frag_peak_end_for_sample = fragment-based peak end coordinate for this sample_name  
frag_peak_width_for_sample = fragment-based peak width (end – start) for this sample_name  

srt_bc_peak_start_for_sample = SRT-barcode-based peak start coordinate for this sample_name  
srt_bc_peak_end_for_sample = SRT-barcode-based peak end coordinate for this sample_name  
srt_bc_peak_width_for_sample = SRT-barcode-based peak width (end – start) for this sample_name  

sample_total_reads_in_consensus_peak = total read count in consensus peak for this sample_name  
sample_total_fragments_in_consensus_peak = total fragment count (merged qbed rows after SRT-barcode correction and deduplication) in consensus peak for this sample_name  
sample_total_insertions_in_consensus_peak = total unique insertions (unique SRT barcodes) in consensus peak for this sample_name  

group_total_reads_in_consensus_peak = sum of sample_total_reads_in_consensus_peak across all sample_name values in this group_name  
group_total_fragments_in_consensus_peak = sum of sample_total_fragments_in_consensus_peak across all sample_name values in this group_name  
group_total_insertions_in_consensus_peak = sum of sample_total_insertions_in_consensus_peak across all sample_name values in this group_name  








