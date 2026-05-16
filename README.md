# Transcriptomic characterization of differentiation dynamics in hiPSC-derived cortical neurons
## Master's Thesis in Life Science Informatics

The thesis addresses cell line imbalance and alternate cell fate commitment in neurons differentiating from human induced pluripotent stem cells by analysing single-cell RNA sequencing data. The codes needed to reproduce the results presented in the thesis are included in this repository.

The original metadata and the variables added in the process of the research are presented in Tables 1 and 2, respectively.

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
| robustID                | Cell name, pool, and sequencing day                              | "AAACCCAAGCAGCCTC_OT_23-1_Pool_d20",                                                                                                                            "AAACCCAAGCGCAATG_OT_23-1_Pool_d20",                                                                                                                            "AAACCCAGTCTTACTT_OT_23-1_Pool_d20"                         |
| predicted.celltype.dfPrimary | Cell type suggested by the first mapping                    | "oRG", "vRG", "vRG"                                         |
| predicted.celltype.score.dfPrimary | Certainty score for the first mapping                 | 0.5255332, 0.5498573, 0.4663895                             |
| mapRef.dfPrimary        | Reference used in the first mapping                              | "primaryMapRef_Polioudakis", "primaryMapRef_Polioudakis",                                                                                                       "primaryMapRef_Polioudakis"                                 |
| predicted.celltype.dfSecondary | Cell type suggested by the second mapping                 | "AstroHindb", "AstroHindb", "panRG"                         |
| predicted.celltype.score.dfSecondary | Certainty score for the second mapping              | 0.5534580, 0.6376691, 0.5051829                             |
| mapRef.dfSecondary      | Reference used in the second mapping                             | "secondaryMapRef_BhaduriOrganoids",                                                                                                                             "secondaryMapRef_BhaduriOrganoids",                                                                                                                             "secondaryMapRef_BhaduriOrganoids"                          |
| twoPassAnnotation_clean | Cell type based on the two-pass annotation                        | "oRG", "vRG", "panRG-O"                             |
| Day_fixed               | Sequencing day rounded to the nearest ten                         | 20, 20, 20                                                  |
| unifiedSampleID         | SampleID in a unified form for all pools                          | "OT_NI23-1_Pool_d20", "OT_NI23-1_Pool_d20", "OT_NI23-1_Pool_d20" |
donorSample             | Columns "SampleID" and "donor" combined                           | "OT_NI23-1_Pool_d20.HPSI0316i-aask_4", "OT_NI23-1_Pool_d20.HPSI0316i-aask_4", "OT_NI23-1_Pool_d20.HPSI0316i-aask_4" |
[^1]: For the HipSci donors, number at the end of the donor name indicated the iPSC line number derived from the donor tissue sample
