# Commands used for Bismark
# (do not run whole script)


# Methylation Extraction (includes bedGraph with scaffolds)
./Bismark-0.22.3/bismark_methylation_extractor -s -o ./IVPF15_MethylExt_bedGraph --bedGraph --scaffolds --ignore 10 IVPF15_ED_1_bismark_bt2_pe.bam


# Only bedGraph (and only on CpG files)
./Bismark-0.22.3/bismark2bedGraph -o IVPF15_bedGraph --dir ./IVPF15_MethylExt_Ign10/bedGraph --scaffolds ./IVPF15_MethylExt_Ign10/CpG*