| Variable name           | Explanation                                                      | First three values                                    |
|-------------------------------------------------|------------------------------------------|-------------------------------------------------------|
| orig.ident              | aa                                                                | "OT_Helsinki", "OT_Helsinki", "OT_Helsinki"                 |
| nCount_RNA              | Total number of molecules in cell                                 | 14613, 22319, 20915                                         |
| nFeature_RNA            | Total number of genes in cell                                     | 4840, 5902, 6120                                            |
| ExperimentID            | Experiment where the data was produced                            | "NI23.1", "NI23.1", "NI23.1"                                |
| SampleID                | Pool name and fixed sequencing day                                | "OT_23-1_Pool_d20", "OT_23-1_Pool_d20", "OT_23-1_Pool_d20"  |
| BatchID                 | Sequencing batch                                                  | "Batch1_090323", "Batch1_090323", "Batch1_090323"           |
| Day                     | Sequencing day                                                    | 20, 20, 20                                                  |
| donor                   | Donor name $^\text{a}$                                             | "HPSI0316i-aask_4", "HPSI0316i-aask_4", "HPSI0516i-oadp_4"  |
| percent.mt              | Percentage of mitochondrial genes                                 | 6.320974, 5.052858, 3.584229                                |
| percent.ribo            | Percentage of ribosomal genes                                     | 13.26447, 18.22254, 16.67861                                |
| Condition               | Patient NDD, KO and correction status                             | "Kabuki", "Kabuki", "Kabuki"                                |
| robustID                | Cell name, pool, and sequencing day                               | "AAACCCAAGCAGCCTC_OT_23-1_Pool_d20",                                                                                                                            "AAACCCAAGCGCAATG_OT_23-1_Pool_d20",                                                                                                                            "AAACCCAGTCTTACTT_OT_23-1_Pool_d20"                         |
| predicted.celltype.dfPrimary | Cell type suggested by the first mapping                     | "oRG", "vRG", "vRG"                                         |
| predicted.celltype.score.dfPrimary | Certainty score for the first mapping                  | 0.5255332, 0.5498573, 0.4663895                             |
| mapRef.dfPrimary        | Reference used in the first mapping                               | "primaryMapRef_Polioudakis", "primaryMapRef_Polioudakis",                                                                                                       "primaryMapRef_Polioudakis"                                 |
| predicted.celltype.dfSecondary | Cell type suggested by the second mapping                  | "AstroHindb", "AstroHindb", "panRG"                         |
| predicted.celltype.score.dfSecondary | Certainty score for the second mapping               | 0.5534580, 0.6376691, 0.5051829                             |
| mapRef.dfSecondary      | Reference used in the second mapping                              | "secondaryMapRef_BhaduriOrganoids",                                                                                                                             "secondaryMapRef_BhaduriOrganoids",                                                                                                                             "secondaryMapRef_BhaduriOrganoids"                          |
| twoPassAnnotation       | Cell type annotated based on the two mappings                     | "oRG", "vRG", "panRG-O"                                     |
| onlyFirstPass           | Result of the fist mapping after accounting for the certainty score | "oRG", "vRG", "Unmapped"                                  |
| onlySecondPass          | Result of the second mapping after accounting for the certainty score | "AstroHindb", "AstroHindb", "panRG"                     |
| twoPassAnnotation_clean | The final two-pass annotation with rare cell types categorized as "Others" | 
| geneMutatedOrCorrected  | The name of the gene with a mutation or a corrected mutation      | "KMT2D", "KMT2D", "KMT2D"                                   | 
| twoPassAnnotation_clean | Cell type based on the two-pass annotation                        | "oRG", "vRG", "panRG-O"                                     |
| Day_fixed               | Sequencing day rounded to the nearest ten                         | 20, 20, 20                                                  |

$^\text{a}$ For the HipSci donors, number at the end of the donor name indicated the iPSC line number derived from the donor tissue sample
