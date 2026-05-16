# Transcriptomic characterization of differentiation dynamics in hiPSC-derived cortical neurons
## Master's Thesis in Life Science Informatics

The thesis addresses cell line imbalance and alternate cell fate commitment in neurons differentiating from human induced pluripotent stem cells by analysing single-cell RNA sequencing data. The codes needed to reproduce the results presented in the thesis are included in this repository.

The original metadata, the additional metadata, and the variables added in the process of the research are presented in Tables 1, 2, and 3, respectively.

**Table 1.** The variables in the original metadata of the single-cell RNA sequencing data object.
| Variable name           | Explanation                                                      | First three values                                    |
|-----------------------------------------|--------------------------------------------------|----------------------------------------|
| orig.ident              | Identity class                                                   | "OT_Helsinki", "OT_Helsinki", "OT_Helsinki"                 |
| nCount_RNA              | Total number of molecules in cell                                | 14613, 22319, 20915                                         |
| nFeature_RNA            | Total number of genes in cell                                    | 4840, 5902, 6120                                            |
| ExperimentID            | Experiment where the data was produced                           | "NI23.1", "NI23.1", "NI23.1"                                |
| SampleID                | Pool name and fixed sequencing day                               | "OT_23-1_Pool_d20", "OT_23-1_Pool_d20", "OT_23-1_Pool_d20"  |
| BatchID                 | Sequencing batch                                                 | "Batch1_090323", "Batch1_090323", "Batch1_090323"           |
| Day                     | Sequencing day                                                   | 20, 20, 20                                                  |
| donor                   | Donor name[^1]                                           | "HPSI0316i-aask_4", "HPSI0316i-aask_4", "HPSI0516i-oadp_4"  |
| percent.mt              | Percentage of mitochondrial genes                                | 6.320974, 5.052858, 3.584229                                |
| percent.ribo            | Percentage of ribosomal genes                                    | 13.26447, 18.22254, 16.67861                                |
| Condition               | Patient NDD, KO and correction status                            | "Kabuki", "Kabuki", "Kabuki"                                |
| robustID                | Cell name, pool, and sequencing day                              | "AAACCCAAGCAGCCTC_OT_23-1_Pool_d20", "AAACCCAAGCGCAATG_OT_23-1_Pool_d20","AAACCCAGTCTTACTT_OT_23-1_Pool_d20" |
| predicted.celltype.dfPrimary | Cell type suggested by the first mapping                    | "oRG", "vRG", "vRG"                                         |
| predicted.celltype.score.dfPrimary | Certainty score for the first mapping                 | 0.5255332, 0.5498573, 0.4663895                             |
| mapRef.dfPrimary        | Reference used in the first mapping                              | "primaryMapRef_Polioudakis", "primaryMapRef_Polioudakis",                                                                                                       "primaryMapRef_Polioudakis"                                 |
| predicted.celltype.dfSecondary | Cell type suggested by the second mapping                 | "AstroHindb", "AstroHindb", "panRG"                         |
| predicted.celltype.score.dfSecondary | Certainty score for the second mapping              | 0.5534580, 0.6376691, 0.5051829                             |
| mapRef.dfSecondary      | Reference used in the second mapping                             | "secondaryMapRef_BhaduriOrganoids", "secondaryMapRef_BhaduriOrganoids", "secondaryMapRef_BhaduriOrganoids"                          |
| twoPassAnnotation_clean | Cell type based on the two-pass annotation                        | "oRG", "vRG", "panRG-O"                             |
| Day_fixed               | Sequencing day rounded to the nearest ten                         | 20, 20, 20                                                  |
| unifiedSampleID         | SampleID in a unified form for all pools                          | "OT_NI23-1_Pool_d20", "OT_NI23-1_Pool_d20", "OT_NI23-1_Pool_d20" |
donorSample             | Columns "SampleID" and "donor" combined                           | "OT_NI23-1_Pool_d20.HPSI0316i-aask_4", "OT_NI23-1_Pool_d20.HPSI0316i-aask_4", "OT_NI23-1_Pool_d20.HPSI0316i-aask_4" |
[^1]: For the HipSci donors, number at the end of the donor name indicated the iPSC line number derived from the donor tissue sample

**Table 2.** The additional, donor-level metadata variables relevant to the analyses.
| Variable name | Explanation | First three values |
|---------------|-------------|-------------------|
| gender | Sex of the donor | "female", "female", "male" |
| disease | Donor NDD status | "Kabuki syndrome", "Kabuki syndrome", "Kabuki syndrome" |
| rnaseq.sendai_reads | Number of reads positive for Sendai virus | 108102, 108102, 6104 |
| pluri_raw | Pluripotency score: how pluripotent the cells are? | 25.555, 25.555, 34.876 | 
| pluri_novelty | Novelty score: how different is the line from normal pluripotent samples? | 1.658, 1.658, 1.408 |


**Table 3.** New variables created and used in the analyses.
| Variable name | Explanation | First three values |
|---------------|-------------|-------------------|
| pool | SampleID with more cohesive names and time point removed | "OT_23-1_Pool", "OT_23-1_Pool", "OT_23-1_Pool" |
| freq_donor | Number of donor cells per cell line | 1825, 1825, 1057 | 
| first_day | The first time point available for the pool | 20, 20, 20 |
| donors_per_pool | The original number of donors in pool | 6, 6, 6 |
| proportion_d0 | Donor cells' proportion of all pool cells on day 0[^1] | 0.1666667, 0.1666667, 0.1666667 |
| proportion_day_fixed | Donor cells' proportion of all pool cells on the sequencing day | 0.2248645, 0.2248645, 0.1302366 |
| R[^2] | R number; reflects cell line imbalance | 1.3491868, 1.3491868, 0.7814194 |
| donor_short | Short form of the donor name | "aask_4", "aask_4", "oadp_4" |
| prop_patient_orig | Proportion of patients in all donors of the pool on day 0[^1] | XXX |
| patient_prop | Proportion of patients in all donors of the pool on the sequencing day |XXX|
| sendai | Is the number of reads positive for sendai virus > 0? | 1, 1, 1 |
| sendai_approx | Is the number of reads positive for sendai virus $\ge$ 10? | 1, 1, 1 |
| proportion_Per[^2] | Cell line pericyte proportion | 0.1063013699, 0.1063013699, 0.0009460738 |
| proportion_vRG[^2] | Cell line vRG cell proportion | 0.4400000, 0.4400000, 0.5600757	 |
| proportion_panRG_O[^2] | Cell line panRG cell proportion | 0.08109589, 0.08109589, 0.14191107 |
| proportion_PgS[^2] | Cell line PgS cell proportion | 0.03452055, 0.03452055, 0.07000946 |
| proportion_oRG[^2] | Cell line oRG cell proportion | 0.10575342, 0.10575342, 0.10575342	|
| proportion_PgG2M[^2] | Cell line PgG2M cell proportion | 0.06684932, 0.06684932, 0.09933775 |
| diff_method | Were the cells differentiated individually, in a cell village, or pooled post-mitotically? | "pool", "pool", "pool" |
| new_SampleID  | SampleID with underscores replaced with hyphens for smoother join with the pseudobulk object | "OT-23-1-Pool-d20", "OT-23-1-Pool-d20", "OT-23-1-Pool-d20"|
| new_donor | donor with underscores replaced with hyphens for smoother join with the pseudobulk object | "HPSI0316i-aask-4", "HPSI0316i-aask-4", "HPSI0516i-oadp-4" |
| cell_type_group | Is the cell a pericyte or of any other cell type? | "nonPericyte", "nonPericyte", "nonPericyte" |
| s[^2][^3] | Cell line level cell cycle score for the synthesis phase | XXX |
| g2m[^2][^3] | Cell line level cell cycle score for the G2 to mitosis transition | XXX |
[^1]: For post-mitotic pools, day 0 is defined as the pooling day.
[^2]: For the indicated variables, the standard scores were also computed and appointed to a variable name consisting of the prefix "z_" and the original variable name.
[^3]: Only included in the pseudobulk objects
