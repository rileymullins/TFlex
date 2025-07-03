# multiplexed_srt_tf_mapping

Identify peaks from multiplexed self-reporting transposon data for identifying transcription factor binding sites.

---

## Annotation file format

This file should be a `.csv` or `.tsv` (tab-separated values) format. The columns must appear **in this order**:  
**library_name**, **sample_barcode**, **sample_name**, **group_name**.

**Example annotation file (tab-separated):**

```
library_name	sample_barcode	sample_name	group_name
Library_A	AAGGCAGACG	Rep1_HyPBase	HyPBase
Library_A	TACCGCTGAC	Rep1_TF_A	TF_A
Library_A	AAGATTAGAC	Rep1_TF_B	TF_B
Library_A	AGTATGACCG	Rep1_TF_C	TF_C
Library_A	TACTAGGTTG	Rep1_TF_D	TF_D
Library_A	CTCTGGATTA	Rep1_TF_E	TF_E
Library_A	AGTTAACGAT	Rep1_TF_F	TF_F
Library_A	GAGACACGTC	Rep1_TF_G	TF_G
Library_A	GTAATGCTTC	Rep1_TF_H	TF_H
Library_A	TTAACGATCG	Rep1_TF_I	TF_I
Library_B	AAGGCAGACG	Rep2_HyPBase	HyPBase
Library_B	TACCGCTGAC	Rep2_TF_A	TF_A
Library_B	AAGATTAGAC	Rep2_TF_B	TF_B
Library_B	AGTATGACCG	Rep2_TF_C	TF_C
Library_B	TACTAGGTTG	Rep2_TF_D	TF_D
Library_B	CTCTGGATTA	Rep2_TF_E	TF_E
Library_B	AGTTAACGAT	Rep2_TF_F	TF_F
Library_B	GAGACACGTC	Rep2_TF_G	TF_G
Library_B	GTAATGCTTC	Rep2_TF_H	TF_H
Library_B	TTAACGATCG	Rep2_TF_I	TF_I
```

- **All columns must be present, and tab-separated.**
- **All `sample_name` values must be unique.**

---

## Variable names

| Variable | Description |
|---|---|
| consensus_peak_id | chrom:consensus_peak_start-consensus_peak_end |
| chrom | chromosome |
| consensus_peak_start | start coordinate of final merged, pan-dataset consensus peak |
| consensus_peak_end | end coordinate of final merged, pan-dataset consensus peak |
| consensus_peak_width | width (end – start) of final merged, pan-dataset consensus peak |
| num_samples_in_consensus_peak | number of unique sample_name values in this consensus peak |
| num_groups_in_consensus_peak | number of unique group_name values in this consensus peak |
| sample_name | sample identifier (from annotation file). **Must be unique.** This is the peak calling unit; could be a replicate (e.g., replicate-1_TF-A), a unique time point (e.g., timepoint1_replicate1_TF-A), or a unique experimental condition (e.g., drugX_replicate1_TF-A). |
| group_name | group identifier (from annotation file). Used to aggregate stats across samples belonging to the same broader group. Examples: "TF_A", "timepoint1_TF-A", or "drugX_TF-A". |
| fragment_peak_start_for_sample | fragment-based peak start coordinate for this sample_name |
| fragment_peak_end_for_sample | fragment-based peak end coordinate for this sample_name |
| fragment_peak_width_for_sample | fragment-based peak width (end – start) for this sample_name |
| srt_bc_peak_start_for_sample | SRT-barcode-based peak start coordinate for this sample_name |
| srt_bc_peak_end_for_sample | SRT-barcode-based peak end coordinate for this sample_name |
| srt_bc_peak_width_for_sample | SRT-barcode-based peak width (end – start) for this sample_name |
| sample_total_reads_in_consensus_peak | total read count in consensus peak for this sample_name |
| sample_total_fragments_in_consensus_peak | total fragment count (merged qbed rows after SRT-barcode correction and deduplication) in consensus peak for this sample_name |
| sample_total_insertions_in_consensus_peak | total unique insertions (unique SRT barcodes) in consensus peak for this sample_name |
| group_total_reads_in_consensus_peak | sum of sample_total_reads_in_consensus_peak across all sample_name values in this group_name |
| group_total_fragments_in_consensus_peak | sum of sample_total_fragments_in_consensus_peak across all sample_name values in this group_name |
| group_total_insertions_in_consensus_peak | sum of sample_total_insertions_in_consensus_peak across all sample_name values in this group_name |

---

