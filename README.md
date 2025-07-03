![Picture1](https://github.com/user-attachments/assets/1966b762-87d9-4af6-91ea-a5500bcd8497)# multiplexed_srt_tf_mapping

Identify peaks from multiplexed self-reporting transposon data for identifying transcription factor binding sites.

---

## Overview of method
This method is capable of antibody-independent semi-multiplxed TF binding site experiments. That is, after initial transfection of a TF-tranposase fusion and a pool of self-reporting transposons containing many SRT barcodes and a unique sample barcode, all samples can be combined. 

This method leverages the self-reporting transposon technology. It is "self-reporting" because it contains a promoter that generates transcripts containing the junction of the terminal repeat and downstream genomic DNA. Thus, the sample barcode and SRT barcode are encoded into the genomic DNA in the transposon, which is then transcribed. Thus, the approximate TF binding site can be identifed through RNA sequencing.

After pooling all samples, the workflow is RNA extraction, cDNA generation, amplification of the self-reporting transpsons, and then library preparation with the Illumina Nextera XT kit.

![image](https://github.com/user-attachments/assets/f67b7cb9-9c70-4061-b525-833f69a9db96)


## Description of data processing script
This pipeline processes multiplexed self-reporting transposon data to identify transcription factor (TF) binding sites. It supports high-throughput experiments with many barcoded samples, corrects technical errors in barcode assignment, collapses redundant reads, calls peaks per sample, and generates a consensus set of reproducible peaks across all experimental conditions.


![Upload<svg width="4913" height="2249" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xml:space="preserve" overflow="hidden"><g transform="translate(32 4)"><g><rect x="18.5001" y="0.499836" width="4834" height="2240" stroke="#000000" stroke-width="4.58333" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="#FFFFFF" fill-opacity="1"/><rect x="1578.5" y="1221.5" width="1160" height="85.0002" stroke="#000000" stroke-width="6.875" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="#DCEAF7" fill-opacity="1"/><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 1778.77 1281)">Sample_A SRT barcode peak</text><path d="M1655.5 1128.5 1655.5 1213.24" stroke="#FF0000" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M2668.5 1128.5 2668.5 1213.24" stroke="#00B0F0" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 1210.78 141)">Load input </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 1460 141)">files</text><path d="M1690.72 175 1690.72 153.25 1629 153.25 1629 109.75 1690.72 109.75 1690.72 88.0002 1766 131.5Z" fill="#000000" fill-rule="evenodd" fill-opacity="1"/><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 1862.27 105)">Correct sample barcodes </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 1941.91 165)">sequencing errors</text><path d="M2610.72 175 2610.72 153.25 2550 153.25 2550 109.75 2610.72 109.75 2610.72 88.0002 2686 131.5Z" fill="#000000" fill-rule="evenodd" fill-opacity="1"/><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2736.08 105)">Generate .parquet files for each </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2925.45 165)">unique sample</text><text fill="#FB0101" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2934.79 481)">Red</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3057.97 481)">=</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3100.93 481)">+ </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3143.9 481)">strand fragments of a unique SRT barcode.</text><text fill="#15B0F1" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2934.79 541)">Blue</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3057.97 541)">= </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3100.93 541)">–</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3142.76 541)">strand fragments of a second unique SRT barcode.</text><text fill="#FB0101" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2934.79 792)">+ strand in the input file means the read moves </text><text fill="#FB0101" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 4073.18 792)">5’</text><text fill="#FB0101" fill-opacity="1" font-family="Wingdings,Wingdings_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 4115 792)"></text><text fill="#FB0101" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 4165.42 792)">3’.</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3099.79 852)">•</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3161.67 852)">The </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3262.5 852)">minimum</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3499.69 852)">position of the SRT barcode across all its fragments in the </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3161.67 911)">peak is the closest to the transposon.</text><text fill="#15B0F1" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2934.79 1030)">•</text><text fill="#15B0F1" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2996.67 1030)">–</text><text fill="#15B0F1" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3038.49 1030)">strand in the input file means reads move </text><text fill="#15B0F1" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 4044.53 1030)">5’</text><text fill="#15B0F1" fill-opacity="1" font-family="Wingdings,Wingdings_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 4086.36 1030)"></text><text fill="#15B0F1" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 4136.77 1030)">3’.</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3099.79 1090)">•</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3161.67 1090)">The </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3262.5 1090)">maximum</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3511.15 1090)">position of the SRT barcode across all its fragments in the </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3161.67 1150)">peak is the closest to the transposon.</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2934.79 1269)">•</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2996.67 1269)">Therefore, SRT barcode peak start and end are defined as follows:</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3099.79 1328)">•</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3161.67 1328)">Peak start </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3416.04 1328)">= </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="italic" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3459.01 1328)">minimum</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3677.29 1328)">of the minimum(+ strand fragment coordinate) and </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3161.67 1388)">maximum(</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3396.57 1388)">-</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3426.93 1388)">strand fragment coordinate)</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3099.79 1448)">•</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3161.67 1448)">Peak end </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3397.14 1448)">= </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="italic" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3440.11 1448)">maximum</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3672.14 1448)">of the minimum(+strand fragment coordinate) and </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3161.67 1507)">maximum(</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3396.57 1507)">-</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3413.18 1507)">strand fragment coordinate)</text><rect x="1452.5" y="520.5" width="1412" height="85.9999" stroke="#000000" stroke-width="6.875" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="#DCEAF7" fill-opacity="1"/><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 1749.84 580)">Sample_A fragment peak (5’</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2483.43 580)">-</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2501.76 580)">3’)</text><path d="M2398.5 427.5 2398.5 512.238" stroke="#00B0F0" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M1474.5 427.5 1474.5 512.238" stroke="#00B0F0" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M2286.5 427.5 2286.5 512.238" stroke="#00B0F0" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M1736.5 427.5 1736.5 512.238" stroke="#FF0000" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M1886.5 427.5 1886.5 512.238" stroke="#FF0000" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M2717.5 427.5 2717.5 512.238" stroke="#FF0000" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M2466.5 427.5 2466.5 512.238" stroke="#FF0000" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M1819.5 427.5 1819.5 512.238" stroke="#00B0F0" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M2806.5 427.5 2806.5 512.238" stroke="#FF0000" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M1628.5 427.5 1628.5 512.238" stroke="#FF0000" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M2641.5 427.5 2641.5 512.238" stroke="#00B0F0" stroke-width="20.625" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><rect x="1578.5" y="2082.5" width="1554" height="86" stroke="#042433" stroke-width="6.875" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="#FFFFB5" fill-opacity="1"/><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2017.33 2142)">Pan</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2118.17 2142)">-</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2136.5 2142)">data consensus peak</text><rect x="1578.5" y="1799.5" width="1160" height="84.9999" stroke="#000000" stroke-width="6.875" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="#DCEAF7" fill-opacity="1"/><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 1778.77 1859)">Sample_A</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2057.46 1859)">SRT barcode peak</text><rect x="2210.5" y="1941.5" width="922" height="84.9999" stroke="#000000" stroke-width="6.875" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="#F2CFEE" fill-opacity="1"/><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2290.6 2000)">Sample_B</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 2571.33 2000)">SRT barcode peak</text><rect x="18.5001" y="0.499836" width="1161" height="2240" stroke="#000000" stroke-width="6.875" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="#F2F2F2" fill-opacity="1"/><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 453.746 484)">Call fragment</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 804.944 484)">-</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 823.277 484)">based peaks </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 734.475 544)">for each sample</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 126.096 1092)">Refine fragment</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 544.325 1092)">-</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 562.658 1092)">based peaks using the </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 94.5856 1152)">SRT barcode position after SRT barcode </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 419.429 1211)">sequencing error correction</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 185.106 1330)">(with 100bp extension on both sides)</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 117.789 1963)">Generate pan</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 468.987 1963)">-</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 487.32 1963)">dataset consensus peaks </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 201.434 2023)">with per</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 412.841 2023)">-</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 431.174 2023)">group and per</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 801.851 2023)">-</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 820.184 2023)">sample stats</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 752.781 133)">Pre</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 841.011 133)">-</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="700" font-stretch="normal" font-size="55" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 859.344 133)">processing</text><path d="M0 0 0.000360892 2239.51" stroke="#000000" stroke-width="6.875" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd" transform="matrix(1 0 0 -1 1179.5 2240.01)"/><path d="M18.5001 274.5 4852.45 274.501" stroke="#000000" stroke-width="6.875" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M18.5001 735.5 4852.45 735.5" stroke="#000000" stroke-width="6.875" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><path d="M18.5001 1720.5 4852.45 1720.5" stroke="#000000" stroke-width="6.875" stroke-linecap="butt" stroke-linejoin="miter" stroke-miterlimit="8" stroke-opacity="1" fill="none" fill-rule="evenodd"/><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3228.24 1790)">Overlapping peaks are merged to allow for downstream differential peak </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3228.24 1850)">analyses and prevent false</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3827.52 1850)">-</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3844.13 1850)">positive conclusions of differential binding of </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3228.24 1909)">slightly offset peaks.</text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3228.24 2028)">The start is defined as the minimum start position of all overlapping </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3228.24 2088)">peaks, and the end is defined as the maximum end position all </text><text fill="#000000" fill-opacity="1" font-family="Arial,Arial_MSFontService,sans-serif" font-style="normal" font-variant="normal" font-weight="400" font-stretch="normal" font-size="50" text-anchor="start" direction="ltr" writing-mode="lr-tb" unicode-bidi="normal" text-decoration="none" transform="matrix(1 0 0 1 3228.24 2148)">overlapping peaks.</text></g></g></svg>ing Picture1.svg…]()





### Features

- **Flexible input:** Accepts multiple `.qbed`, `.bed`, or gzipped equivalents from a directory.
- **Barcode correction:** Corrects sample barcodes against a whitelist (from the annotation file) using defined Hamming distance cutoff.
- **Customizable:** Uses a user-supplied annotation file to map barcodes to sample and group names.
- **Fragment deduplication:** Collapses reads within each peak that are identical in coordinates and SRT barcodes using UMI-tools clustering-based correction.
- **Parallel peak calling:** Efficiently calls peaks per sample using pycallingcards, leveraging multiprocessing.
- **Consensus peak construction:** Merges sample-specific peaks to produce a unified, non-redundant set of consensus peaks across the dataset.
- **Comprehensive reporting:** Outputs per-sample and per-group statistics, final peak sets, and detailed column definitions.
- **Robust logging:** Logs all steps and errors for traceability.
- **Memory efficient:** Capable of running on most laptops, but will need to monitor memory usage and lower the number of workers as needed.

---

## Usage

```bash
python ccaller_to_assign_reproducible_peaks_with_full_pan_barcode_consensus_peaks_and_umi_based_peak_boundaries.py \
  --input_dir <input_qbed_dir> \
  --output_dir <output_dir> \
  --annotation_file <annotation_file.tsv or .csv> \
  --workers <num>] \
  --sample_barcode_dist_threshold <int> \
  --srt_bc_dist_threshold <int> \
  [other peak calling parameters, see below]
```

### Required arguments

- `--input_dir`  
  Directory containing input `.qbed`, `.bed`, `.qbed.gz`, or `.bed.gz` files. Each file should be in standard qbed/bed format.
- `--output_dir`  
  Directory where output files will be written.
- `--annotation_file`  
  Path to the annotation file (see format below).

### Optional arguments

- `--workers`  
  Number of parallel worker processes for peak calling (default: 10). One worker will be assigned a sample_name and perform all functions on that sample_name.
- `--sample_barcode_dist_threshold`  
  Maximum Hamming distance allowed for correcting sample barcodes (default: 2). A value of 2 means that only sample barcodes that have ≤ 2 mismatches from the sample barcode are assigned to the whitelisted sample barcode.
- `--srt_bc_dist_threshold`  
  Maximum Hamming distance for SRT-barcode clustering (default: 1). A value of 1 means that only SRT-barcodes with ≤ 1 mismatches from each other will be grouped.
- Additional advanced parameters for peak calling (see script source for defaults and descriptions):
  - `--pvalue_cutoff`, `--pvalue_adj_cutoff`, `--min_insertions`, `--minlen`, `--extend`, `--maxbetween`, `--minnum`, `--lam_win_size`
  - It is recommended to leave these parameters at the default values.
  - The purpose of these parameters is to define all regions of concentrated fragments for subsequent SRT-barcode-based peak boundary definitions and donwstream analysis.
  - Fragment peaks with <5 fragments following SRT barcode correction are considered noise and discarded from the analysis.

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

## Input data requirements

### Supported input file formats

- `.qbed`, `.bed`, `.qbed.gz`, `.bed.gz`
- Each file must contain columns in the order:  
  `chrom`, `start`, `end`, `reads`, `strand`, `name`
- The `name` field should be formatted as:  
  `library_name/sample_barcode/srt_barcode`
- The sample_barcd
  
### Example qbed row

```
chr1	12345	12350	100	+	[library_name]/[sample_barcode]/[srt_barcode]
```

### Interpretation of 'strand', 'start', and 'end' values 

- Note that by convention of the alignment pipeline:
  - The **'+'** strand **'end' coordinate** is the true end of the read.
  - The **'-'** strand **'start' coordinate** is the true end of the read.
    - **'End of the read'** means last base read from the sequencer.
  
  - The strand in the qbed refers the strand that the **transposon** was inserted into, not the strand of the read.
    - This means that the **R2 (i7) read** is for **'+'** strand qbed rows are moving 5' <- 3'
    - The **R2 (i7) read** for **'-'** qbed rows are moving 5' -> 3'.
    - Yes, this is opposite of the conventional pairing of a strand and direction of the read because the strand in the qbed refers to strand that the **transposon** was inserted into, not the strand of the read.
    
- In summary:
  - For a **'+'** strand row in the qbed, the **transposon** is inserted on the **'+'** strand, the **R2 read** is moving  5' <- 3', and the **'end' coordinate** is the true end of the read.
  - For a **'-'** strand row in the qbed, the **transposon** is inserted on the **'-'** strand, the **R2 read** is moving  5' -> 3', and the **'start' coordinate** is the true end of the read.


---

## Pipeline steps

1. **Input loading and barcode parsing:**  
   Reads all qbed/bed files and parses the `name` field into `library_name`, `sample_barcode_raw`, and `srt_bc`.
2. **Barcode correction:**  
   Corrects sample barcodes against the annotation whitelist, using Hamming distance up to the `--sample_barcode_dist_threshold`.
3. **Annotation:**  
   Joins corrected reads with the annotation file, assigning `sample_name` and `group_name`.
4. **Per-sample partitioning:**  
   Collapses all reads for each sample into a `.parquet` file for efficient downstream processing.
5. **Peak calling:**  
   For each sample, calls peaks with pycallingcards, then refines peak boundaries using strand-aware logic of the position per SRT barcode that is the most proximal to the junction of the transposon and genomic DNA.
6. **Consensus peak generation:**  
   Merges sample-specific peaks to produce consensus peaks across all samples.
7. **Read-to-peak mapping:**  
   Intersects all reads and all sample peaks with consensus peaks to generate per-sample and per-group statistics.
8. **Output aggregation:**  
   Produces a final report (`final_results.tsv`) with detailed statistics and peak boundaries, and a `column_definitions.tsv` file describing all output columns.

---

## Output files

- `final_results.tsv`  
  Tab-separated table containing all peak and sample statistics. See variable definitions below.
- `column_definitions.tsv`  
  Tab-separated table describing all columns in `final_results.tsv`.
- `pipeline.log`  
  Log file of all steps, warnings, and errors.
- `collapsed_per_sample/*.parquet`  
  Per-sample data files (intermediate).
- `pybedtools_tmp/`  
  Temporary files for bedtools operations (auto-cleaned at completion).

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

## Dependencies

- Python 3.8+
- [polars](https://www.pola.rs/) (for fast DataFrame operations)
- [pandas](https://pandas.pydata.org/)
- [pybedtools](https://daler.github.io/pybedtools/)
- [pycallingcards](https://github.com/brwnj/pycallingcards)
- [UMI-tools](https://umi-tools.readthedocs.io/en/latest/)
- [tqdm](https://tqdm.github.io/) (for progress bars)

Install dependencies with pip:

```bash
pip install polars pandas pybedtools pycallingcards umi_tools tqdm
```

> Note: pybedtools requires bedtools (>=2.27.1) installed on your system and available in your `$PATH`.

---

## Logging and troubleshooting

- Progress, errors, and warnings are logged to `pipeline.log` in the output directory.
- Intermediate files are saved in the output directory for reproducibility and debugging.
- The script exits with informative error messages if any required step fails.

---

## Advanced options

- All pycallingcards peak calling parameters can be overridden via the command line.
- The pipeline can be parallelized across as many CPUs as available using `--workers`.
- Barcode and SRT clustering stringency can be tuned via `--sample_barcode_dist_threshold` and `--srt_bc_dist_threshold`.

---

## Example run

```bash
python ccaller_to_assign_reproducible_peaks_with_full_pan_barcode_consensus_peaks_and_umi_based_peak_boundaries.py \
  --input_dir /path/to/qbed_files \
  --output_dir /path/to/results \
  --annotation_file /path/to/annotation_file.tsv \
  --workers 8 \
  --sample_barcode_dist_threshold 2 \
  --srt_bc_dist_threshold 1
```

---

## Contact

For questions, bugs, or feature requests, please open an issue or contact the repository maintainer.
