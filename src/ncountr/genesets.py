"""Gene set collections, GMT file parsing, and Hallmark gene sets.

Includes a curated subset of MSigDB Hallmark gene sets (Liberzon et al.
2015) most relevant to nCounter immune and host-response panels.  Gene
symbols sourced from MSigDB v2024.1 (CC BY 4.0).
"""

from __future__ import annotations

from pathlib import Path
from typing import Union


# ---------------------------------------------------------------------------
# MSigDB Hallmark gene sets (curated subset for nCounter panels)
# ---------------------------------------------------------------------------
# Each list contains the core member genes.  The filter_gene_sets()
# function intersects these with the genes actually measured on the panel.

_HALLMARK_SETS: dict[str, list[str]] = {
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": [
        "ADAR", "BST2", "CMPK2", "CNP", "DDX60", "DHX58", "EIF2AK2",
        "EPSTI1", "GBP4", "HERC5", "HERC6", "HLA-C", "IFI27", "IFI30",
        "IFI35", "IFI44", "IFI44L", "IFI6", "IFIH1", "IFIT1", "IFIT2",
        "IFIT3", "IFIT5", "IFITM1", "IFITM2", "IFITM3", "IRF1", "IRF2",
        "IRF7", "IRF9", "ISG15", "ISG20", "LAMP3", "LAP3", "LGALS3BP",
        "LY6E", "MOV10", "MX1", "MX2", "NMI", "NUB1", "OAS1", "OAS2",
        "OAS3", "OASL", "OGFR", "PARP12", "PARP14", "PARP9", "PLSCR1",
        "PNPT1", "PROCR", "PSMA2", "PSMA3", "PSMB8", "PSMB9", "PSME1",
        "PSME2", "RIPK2", "RNF31", "RSAD2", "RTP4", "SAMD9", "SAMD9L",
        "SP110", "STAT1", "STAT2", "TDRD7", "TMEM140", "TRAFD1", "TRIM14",
        "TRIM21", "TRIM25", "TRIM5", "TXNIP", "UBA7", "UBE2L6", "USP18",
        "XAF1",
    ],
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": [
        "ADAR", "APOL6", "B2M", "BATF2", "BST2", "BTG1", "CASP1",
        "CASP3", "CASP4", "CASP7", "CASP8", "CCL2", "CCL5", "CCL7",
        "CD274", "CD38", "CD40", "CD74", "CDKN1A", "CIITA", "CMKLR1",
        "CMPK2", "CSF2RB", "CXCL10", "CXCL11", "CXCL9", "DDX58",
        "DDX60", "DHX58", "EIF2AK2", "EPSTI1", "FAS", "FCGR1A",
        "FGL2", "GBP1", "GBP2", "GBP4", "GBP5", "HLA-A", "HLA-B",
        "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1",
        "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-E",
        "HLA-F", "HLA-G", "ICAM1", "IDO1", "IFI27", "IFI30", "IFI35",
        "IFI44", "IFI44L", "IFIH1", "IFIT1", "IFIT2", "IFIT3",
        "IFITM1", "IFITM2", "IFITM3", "IL10RA", "IL15", "IL15RA",
        "IL18BP", "IL2RB", "IL4R", "IL6", "IL7", "IRF1", "IRF2",
        "IRF4", "IRF5", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20",
        "JAK2", "LAP3", "LGALS3BP", "LY6E", "LYSE", "MARCHF1", "MT2A",
        "MX1", "MX2", "MYD88", "NCOA3", "NFKB1", "NFKBIA", "NMI",
        "NOD1", "NUB1", "OAS1", "OAS2", "OAS3", "OASL", "PARP12",
        "PARP14", "PARP9", "PDE4B", "PELI1", "PML", "PNP", "PNPT1",
        "PSMA2", "PSMA3", "PSMB10", "PSMB2", "PSMB8", "PSMB9",
        "PSME1", "PSME2", "PTGS2", "PTPN1", "PTPN2", "PTPN6",
        "RAPGEF6", "RBCK1", "RIPK2", "RNF31", "RSAD2", "RTP4",
        "SAMD9L", "SECTM1", "SELP", "SERPING1", "SOCS1", "SOCS3",
        "SOD2", "SP110", "SPPL2A", "SRI", "STAT1", "STAT2", "STAT3",
        "STAT4", "TAP1", "TAP2", "TAPBP", "TDRD7", "TNFAIP2",
        "TNFAIP3", "TNFAIP6", "TNFSF10", "TRAFD1", "TRIM14", "TRIM21",
        "TRIM25", "TRIM26", "TXNIP", "UBE2L6", "USP18", "VCAM1",
        "VAMP5", "VAMP8", "XAF1", "XCL1", "ZBP1",
    ],
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": [
        "ABCA1", "ACKR3", "AREG", "ATF3", "B4GALT1", "B4GALT5",
        "BCL2A1", "BCL3", "BCL6", "BIRC2", "BIRC3", "BMP2", "BTG2",
        "CCL2", "CCL20", "CCL4", "CCL5", "CCN1", "CCND1", "CCNL1",
        "CCRL2", "CD44", "CD69", "CD80", "CD83", "CDKN1A", "CEBPB",
        "CEBPD", "CFLAR", "CLCF1", "CSF1", "CSF2", "CXCL1", "CXCL10",
        "CXCL11", "CXCL2", "CXCL3", "CXCL6", "DUSP1", "DUSP2",
        "DUSP4", "DUSP5", "EDN1", "EFNA1", "EGR1", "EGR2", "EGR3",
        "ETS2", "F3", "FJX1", "FOS", "FOSB", "FOSL1", "FOSL2", "FUT4",
        "G0S2", "GADD45A", "GADD45B", "GCH1", "GEM", "GFPT2", "GPR183",
        "HBEGF", "ICAM1", "ICOSLG", "ID2", "IER2", "IER3", "IER5",
        "IFIH1", "IFNGR2", "IKBKE", "IL12B", "IL15RA", "IL18", "IL1A",
        "IL1B", "IL23A", "IL6", "IL7R", "INHBA", "IRF1", "IRS2",
        "JAG1", "JUN", "JUNB", "KDM6B", "KLF10", "KLF2", "KLF4",
        "KLF6", "KLF9", "KYNU", "LAMB3", "LIF", "LITAF", "MAFF",
        "MAP2K3", "MAP3K8", "MARCKS", "MCL1", "MSC", "MXD1", "MYC",
        "NFAT5", "NFE2L2", "NFIL3", "NFKB1", "NFKB2", "NFKBIA",
        "NFKBIE", "NINJ1", "NR4A1", "NR4A2", "NR4A3", "OLR1",
        "PANX1", "PDE4B", "PDLIM5", "PFKFB3", "PHLDA1", "PHLDA2",
        "PLAU", "PLAUR", "PLEK", "PLK2", "PLPP3", "PMEPA1", "PNRC1",
        "PPP1R15A", "PTGER4", "PTGS2", "PTX3", "RCAN1", "REL", "RELA",
        "RELB", "RHOB", "RIPK2", "RNF19B", "SAT1", "SDC4", "SERPINB2",
        "SERPINB8", "SERPINE1", "SGK1", "SIK1", "SLC16A6", "SLC2A3",
        "SLC2A6", "SMAD3", "SNN", "SOCS3", "SOD2", "SPHK1", "SQSTM1",
        "STAR", "TANK", "TAP1", "TGIF1", "TIPARP", "TLR2", "TNAIP3",
        "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFAIP8", "TNFRSF9",
        "TNFSF9", "TNIP1", "TNIP2", "TRAF1", "TRIB1", "TRIP10",
        "TSC22D1", "TUBB2A", "VEGFA", "ZBTB10", "ZC3H12A", "ZFP36",
    ],
    "HALLMARK_IL6_JAK_STAT3_SIGNALING": [
        "A2M", "ACVR1B", "AHR", "ANGPT1", "BAK1", "BCL2L1", "CBL",
        "CCL7", "CCR1", "CD14", "CD36", "CD44", "CD9", "CNTFR", "CSF1",
        "CSF2RB", "CSF3R", "CXCL1", "DNTT", "FAS", "GRB2", "HAX1",
        "HMOX1", "IFNAR1", "IFNGR1", "IFNGR2", "IL10", "IL10RB",
        "IL13RA1", "IL15RA", "IL17RA", "IL1R1", "IL2RA", "IL2RG",
        "IL3RA", "IL4R", "IL6", "IL6R", "IL6ST", "IL7", "IL9R",
        "IRF1", "IRF9", "ISG15", "ITGA4", "JAK1", "JAK2", "JAK3",
        "JUN", "JUNB", "LEPR", "LIF", "LIFR", "MAP3K8", "MYC",
        "OSMR", "PF4", "PIK3R5", "PIM1", "PLA2G4A", "PIAS1", "REL",
        "SOCS1", "SOCS3", "STAT1", "STAT2", "STAT3", "TBX21", "TGFB1",
        "TLR2", "TNFRSF12A", "TNFRSF1A", "TNFRSF1B", "TNFRSF21",
        "TRAF5", "TYK2",
    ],
    "HALLMARK_INFLAMMATORY_RESPONSE": [
        "ABI1", "ACVR1B", "ACVR2A", "ADGRE1", "AHR", "BEST1", "BST1",
        "BTG2", "C3AR1", "C5AR1", "CASP1", "CASP4", "CCL11", "CCL13",
        "CCL17", "CCL2", "CCL20", "CCL22", "CCL24", "CCL3", "CCL4",
        "CCL5", "CCL7", "CCR1", "CCR2", "CCR4", "CCR5", "CCR7", "CD14",
        "CD40", "CD44", "CD48", "CD69", "CD74", "CD82", "CDKN1A",
        "CEBPB", "CLEC5A", "CSF1", "CSF2RB", "CSF3R", "CXCL1", "CXCL10",
        "CXCL11", "CXCL2", "CXCL3", "CXCL6", "CXCL8", "CXCL9",
        "CXCR4", "DCBLD2", "EBI3", "EIF2AK2", "EMP3", "ETS1", "ETS2",
        "F3", "FAS", "FCGR2A", "FCGR2B", "FOS", "FPR1", "GCH1",
        "GPR68", "HIF1A", "HRH1", "ICAM1", "ICOSLG", "IDO1", "IFIH1",
        "IFNGR2", "IL10", "IL10RA", "IL15", "IL15RA", "IL18", "IL18R1",
        "IL1A", "IL1B", "IL1R1", "IL1R2", "IL1RN", "IL2RA", "IL4R",
        "IL6", "IL6R", "IL7", "IL7R", "IRF1", "IRF4", "IRF5", "IRF7",
        "IRAK2", "ITGA5", "ITGB3", "ITGB8", "JAG1", "KLF6",
        "KYNU", "LAMP3", "LCP2", "LMNA", "LYN", "MAP2K3", "MAP3K8",
        "MEFV", "MMP14", "MMP9", "MSR1", "MYC", "MYD88", "NAMPT",
        "NFKB1", "NFKB2", "NFKBIA", "NLRP3", "NMI", "NOD2", "NR4A1",
        "OLR1", "OSMR", "P2RY2", "PANX1", "PCDH7", "PDE4B", "PDPN",
        "PIK3R5", "PLAU", "PLAUR", "PLEK", "PNP", "PPBP", "PROK2",
        "PTAFR", "PTGER2", "PTGER4", "PTGS2", "PTPRE", "RELA", "RELB",
        "RIPK2", "RNF19B", "SCAMP5", "SELL", "SERPINE1", "SLAMF1",
        "SLC11A1", "SLC28A2", "SLC7A2", "SOCS3", "SOD2", "SPI1",
        "SPHK1", "SRI", "STAB1", "TLR1", "TLR2", "TLR3", "TLR4",
        "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF10", "TNFSF14",
        "TNFSF15", "TNFSF4", "TRAF1", "VCAM1", "VEGFA", "ZC3H12A",
    ],
    "HALLMARK_COMPLEMENT": [
        "A2M", "ACPP", "C1QA", "C1QB", "C1QC", "C1R", "C1S", "C2",
        "C3", "C3AR1", "C4A", "C4B", "C4BPA", "C5", "C5AR1", "C6",
        "C7", "C8A", "C8B", "C8G", "C9", "CBL", "CD14", "CD163",
        "CD36", "CD46", "CD55", "CD59", "CD63", "CD81", "CD93",
        "CFB", "CFD", "CFH", "CFI", "CLU", "CR1", "CR2", "CTSD",
        "F10", "F12", "F2", "F2R", "F3", "F5", "F7", "F8", "F9",
        "FGA", "FGB", "FGG", "GNAI2", "GNG12", "HMOX1", "ITGAM",
        "ITGAX", "ITGB2", "KNG1", "LMAN1", "LRP1", "MASP1", "MASP2",
        "MBL2", "PIGA", "PLAT", "PLAU", "PLAUR", "PLG", "PROC",
        "PROS1", "SERPINA1", "SERPINA5", "SERPINC1", "SERPIND1",
        "SERPINE1", "SERPINF2", "SERPING1", "TFPI", "TFPI2", "THBD",
        "TIMP1", "VWF",
    ],
    "HALLMARK_IL2_STAT5_SIGNALING": [
        "AHR", "BATF", "BCL2", "BCL2L1", "CASP3", "CCND2", "CCND3",
        "CCR4", "CD25", "CD4", "CD44", "CD69", "CD80", "CD83",
        "CDKN1A", "CISH", "CSF1", "CSF2", "CSF2RB", "CTLA4", "EGR1",
        "ETS1", "FOXP3", "GATA3", "GPI", "GZMB", "HIF1A", "HMOX1",
        "ICOS", "IFNG", "IKZF2", "IL10", "IL10RA", "IL12RB2", "IL13",
        "IL15RA", "IL17A", "IL17F", "IL1R1", "IL1R2", "IL1RL1", "IL2",
        "IL21", "IL21R", "IL2RA", "IL2RB", "IL2RG", "IL3", "IL4",
        "IL4R", "IL5", "IL6R", "IL7R", "IL9", "IL9R", "IRF4", "ITGA6",
        "ITGAE", "JAK1", "JAK3", "JUN", "KLF2", "LMNA", "LTA",
        "MAPKAPK2", "MYC", "NFIL3", "PENK", "PIM1", "PRDM1",
        "PTGER2", "PTGER4", "SELL", "SH2D1A", "SOCS1", "SOCS2",
        "SOCS3", "STAT1", "STAT3", "STAT4", "STAT5A", "STAT5B",
        "TBX21", "TNFRSF18", "TNFRSF4", "TNFRSF9", "TNFSF11",
        "TNFSF4", "TYK2",
    ],
    "HALLMARK_ALLOGRAFT_REJECTION": [
        "AKT1", "B2M", "BATF", "C2", "CAPG", "CASP1", "CASP4",
        "CCL11", "CCL13", "CCL19", "CCL2", "CCL22", "CCL4", "CCL5",
        "CCL7", "CCR1", "CCR2", "CCR5", "CD2", "CD247", "CD28",
        "CD3D", "CD3E", "CD3G", "CD4", "CD40", "CD40LG", "CD47",
        "CD7", "CD74", "CD80", "CD86", "CD8A", "CD8B", "CDKN2A",
        "CSF2", "CTLA4", "CXCL10", "CXCL11", "CXCL13", "CXCL9",
        "CXCR3", "ETS1", "F2R", "FAS", "FASLG", "FCGR2A", "FCGR2B",
        "FYN", "GATA3", "GBP2", "GNLY", "GZMA", "GZMB", "HLA-A",
        "HLA-B", "HLA-C", "HLA-DMA", "HLA-DMB", "HLA-DPA1", "HLA-DPB1",
        "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-E",
        "ICAM1", "ICOS", "IFNG", "IFNGR1", "IKBKB", "IL10", "IL11",
        "IL12A", "IL12B", "IL12RB1", "IL13", "IL15", "IL16", "IL18",
        "IL1B", "IL2", "IL2RA", "IL2RB", "IL2RG", "IL4", "IL4R",
        "IL6", "IL7", "IRF1", "IRF4", "IRF7", "IRF8", "ISG20",
        "ITGAL", "ITGB2", "JAK2", "KLF2", "LAMP3", "LCK", "LCP2",
        "LIF", "LTA", "LTB", "MAP3K7", "MYD88", "NFKB1", "NFKB2",
        "NFKBIA", "NOS2", "OAS1", "PF4", "PIK3CA", "PML", "PRF1",
        "PRKCB", "PSMB9", "PTPN6", "REL", "RELA", "RIPK2", "SEMA4D",
        "SOCS1", "SPI1", "STAT1", "STAT4", "TAP1", "TAP2", "TAPBP",
        "TBX21", "TLR1", "TLR2", "TLR3", "TLR4", "TLR6", "TNF",
        "TNFRSF1A", "TNFRSF1B", "TNFSF10", "VCAM1", "ZAP70",
    ],
    "HALLMARK_APOPTOSIS": [
        "ADD1", "AHR", "ANXA1", "APP", "ATF3", "BAD", "BAK1", "BAX",
        "BCL10", "BCL2", "BCL2L1", "BCL2L11", "BID", "BIRC3", "BMF",
        "BMP2", "BNIP3L", "BRCA1", "BTG2", "BTG3", "CASP1", "CASP3",
        "CASP4", "CASP6", "CASP7", "CASP8", "CASP9", "CCND1", "CCND2",
        "CDKN1A", "CFLAR", "CLU", "CREBBP", "CYCS", "DDIT3", "DFFA",
        "DIABLO", "DNAJA1", "DNAJC3", "EBP", "EGR3", "EREG", "ETF1",
        "F2R", "FADD", "FAS", "FDXR", "FEZ1", "GADD45A", "GADD45B",
        "GCH1", "GNA15", "GPX1", "GPX3", "GPX4", "GSN", "GUCY2D",
        "HERPUD1", "HMGB2", "HSPB1", "IER3", "IFNB1", "IGF2R",
        "IL1A", "IL1B", "IL6", "ISG20", "JUN", "LGALS3", "LMNA",
        "MCL1", "MYC", "NEDD9", "NFKB1", "NFKBIA", "PEA15", "PLAT",
        "PLCB2", "PPT1", "PSEN1", "PTPN13", "RELA", "RETSAT", "RHOB",
        "RNASEL", "ROCK1", "SAT1", "SATB1", "SCN9A", "SMAD7", "SOD1",
        "SOD2", "SPTAN1", "SQSTM1", "TAP1", "TGFB1", "TGFBR3",
        "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF10B", "TNFRSF12A",
        "TNFRSF1A", "TNFSF10", "TOP2A", "TP53", "TSPO",
    ],
    "HALLMARK_COAGULATION": [
        "A2M", "ACOX2", "ADAM9", "C1QA", "C1R", "C1S", "C2", "C3",
        "C4BPA", "C8A", "C8G", "C9", "CD36", "CD46", "CD55", "CD59",
        "CD63", "CD9", "CFB", "CFH", "CLU", "CPB2", "CTSE", "F10",
        "F11", "F12", "F13A1", "F2", "F2R", "F3", "F5", "F7", "F8",
        "F9", "FGA", "FGB", "FGG", "FN1", "GDA", "GP1BA", "GP5",
        "GP6", "GP9", "HMOX1", "ITGA2B", "ITGB3", "KNG1", "LMAN1",
        "MASP1", "MASP2", "MBL2", "MMP1", "MMP14", "MMP2", "MMP3",
        "MMP9", "PDGFA", "PF4", "PLA2G4A", "PLAT", "PLAU", "PLAUR",
        "PLG", "PROC", "PROS1", "SERPINA1", "SERPINA10", "SERPINA5",
        "SERPINB2", "SERPINC1", "SERPIND1", "SERPINE1", "SERPINF2",
        "SERPING1", "TFPI", "TFPI2", "THBD", "THBS1", "TIMP1",
        "TIMP3", "VWF",
    ],
    "HALLMARK_HYPOXIA": [
        "ACKR3", "ADM", "ADORA2B", "AK4", "ALDOA", "ALDOB", "ALDOC",
        "ANGPTL4", "ANKZF1", "ANXA2", "ATF3", "BCL2", "BGN", "BNIP3",
        "BNIP3L", "CA9", "CASP6", "CAVIN1", "CCNG2", "CDKN1A",
        "CDKN1B", "CDKN1C", "CITED2", "COL5A1", "CP", "CREBBP",
        "DDIT3", "DDIT4", "DTNA", "EDN1", "EFNA1", "EFNA3", "ENO1",
        "ENO2", "ENO3", "ERO1A", "ETS1", "EXT1", "F3", "FAM162A",
        "FOS", "FOSL2", "GAA", "GAPDH", "GBE1", "GCK", "GLRX",
        "GPC1", "GPI", "GYS1", "HIF1A", "HK1", "HK2", "HMOX1",
        "IER3", "IGFBP1", "IGFBP3", "IL6", "IRS2", "ISG20", "JUN",
        "KDM3A", "KLF6", "KLF7", "LDHA", "LDHB", "LOX", "MXI1",
        "NAGK", "NDRG1", "NFIL3", "NR3C1", "P4HA1", "P4HA2",
        "PCDH7", "PDK1", "PDK3", "PFKFB3", "PFKL", "PFKP", "PGK1",
        "PGM1", "PGM2", "PKLR", "PKM", "PLAUR", "PLOD1", "PLOD2",
        "POU5F1", "PPP1R15A", "PRDX5", "RORA", "SAT1", "SCARB1",
        "SDC2", "SDC3", "SELENBP1", "SERPINE1", "SLC16A3", "SLC2A1",
        "SLC2A3", "SLC6A6", "STC1", "STC2", "TGFB3", "TGFBI",
        "TGM2", "TIPARP", "TKTL1", "TMEM45A", "TNFAIP3", "TPI1",
        "VEGFA", "XPNPEP1",
    ],
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": [
        "ABI3BP", "ACTA2", "ADAM12", "ANPEP", "APLP1", "AREG", "BASP1",
        "BGN", "BMP1", "CAP2", "CAPG", "CD44", "CD59", "CDH11",
        "CDH2", "CDH6", "COL11A1", "COL12A1", "COL16A1", "COL1A1",
        "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2",
        "COL5A3", "COL6A2", "COL6A3", "COL7A1", "COL8A2", "COMP",
        "CRLF1", "CTHRC1", "CXCL1", "CXCL12", "CXCL6", "CXCL8",
        "DKK1", "DSP", "ECM1", "ECM2", "EDIL3", "EFEMP2", "ELN",
        "EMP3", "ENO2", "FAP", "FBLN1", "FBLN2", "FBLN5", "FBN1",
        "FBN2", "FERMT2", "FGF2", "FLNA", "FMOD", "FN1", "FOXC2",
        "FSTL1", "FSTL3", "GAS1", "GEM", "GJA1", "GLIPR1", "GPC1",
        "GPX7", "GREM1", "HTRA1", "ICAM1", "ID2", "IGFBP2", "IGFBP3",
        "IGFBP4", "IL15", "IL32", "IL6", "INHBA", "ITGA2", "ITGA5",
        "ITGAV", "ITGB1", "ITGB3", "ITGB5", "JUN", "LAMA1", "LAMA2",
        "LAMA3", "LAMB3", "LAMC1", "LAMC2", "LGALS1", "LOX", "LOXL1",
        "LOXL2", "LRP1", "LRRC15", "LTBP1", "LTBP2", "MCM7", "MGP",
        "MMP1", "MMP14", "MMP2", "MMP3", "MSX1", "MXRA5", "NOTCH2",
        "NT5E", "NTM", "PCOLCE", "PCOLCE2", "PDGFRB", "PFN2",
        "PLAUR", "PLOD1", "PLOD2", "PLOD3", "PMEPA1", "PMP22",
        "POSTN", "PPIB", "PRRX1", "PTHLH", "PTX3", "PVR", "RGS4",
        "RHOB", "SAT1", "SCG2", "SDC1", "SDC4", "SERPINE1", "SERPINE2",
        "SERPINH1", "SFRP1", "SFRP4", "SLC6A8", "SLIT2", "SLIT3",
        "SNAI2", "SPARC", "SPOCK1", "SPP1", "TAGLN", "TGM2",
        "TGFB1", "TGFBI", "THBS1", "THBS2", "THY1", "TIMP1", "TIMP3",
        "TNC", "TNFAIP3", "TNFRSF11B", "TNFRSF12A", "TPM1", "TPM2",
        "TPM4", "VCAM1", "VEGFA", "VEGFC", "VIM", "WIPF1", "WNT5A",
    ],
}


# ---------------------------------------------------------------------------
# GMT file parsing
# ---------------------------------------------------------------------------

def load_gmt(path: Union[str, Path]) -> dict[str, list[str]]:
    """Load gene sets from a GMT file.

    GMT format: ``set_name<TAB>description<TAB>gene1<TAB>gene2<TAB>...``

    Parameters
    ----------
    path : str or Path
        Path to ``.gmt`` file.

    Returns
    -------
    dict[str, list[str]]
    """
    gene_sets: dict[str, list[str]] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = [g for g in parts[2:] if g]
            gene_sets[name] = genes
    return gene_sets


def save_gmt(
    gene_sets: dict[str, list[str]],
    path: Union[str, Path],
    *,
    descriptions: dict[str, str] | None = None,
) -> None:
    """Write gene sets to GMT format."""
    descs = descriptions or {}
    with open(path, "w") as fh:
        for name, genes in gene_sets.items():
            desc = descs.get(name, "")
            fh.write(f"{name}\t{desc}\t" + "\t".join(genes) + "\n")


# ---------------------------------------------------------------------------
# Lookup functions
# ---------------------------------------------------------------------------

def get_hallmark_set(name: str) -> list[str]:
    """Return a Hallmark gene set by name (case-insensitive).

    Accepts names with or without the ``HALLMARK_`` prefix.
    """
    key = name.upper()
    if not key.startswith("HALLMARK_"):
        key = "HALLMARK_" + key
    for k, v in _HALLMARK_SETS.items():
        if k == key:
            return list(v)
    raise KeyError(
        f"Unknown Hallmark set {name!r}. "
        f"Available: {list(_HALLMARK_SETS.keys())}"
    )


def list_hallmark_sets() -> list[str]:
    """Return names of all built-in Hallmark gene sets."""
    return list(_HALLMARK_SETS.keys())


def get_all_gene_sets() -> dict[str, list[str]]:
    """Return merged dict of all built-in gene sets (custom + Hallmark)."""
    from ncountr.datasets import _GENE_SETS
    merged = dict(_GENE_SETS)
    merged.update(_HALLMARK_SETS)
    return merged


def filter_gene_sets(
    gene_sets: dict[str, list[str]],
    measured_genes: list[str],
    *,
    min_overlap: int = 5,
    max_overlap: int = 500,
) -> dict[str, list[str]]:
    """Filter gene sets to those with sufficient overlap with measured genes.

    Returns a new dict where each gene list contains only the genes
    that appear in *measured_genes*.
    """
    measured = set(measured_genes)
    filtered: dict[str, list[str]] = {}
    for name, genes in gene_sets.items():
        overlap = [g for g in genes if g in measured]
        if min_overlap <= len(overlap) <= max_overlap:
            filtered[name] = overlap
    return filtered
