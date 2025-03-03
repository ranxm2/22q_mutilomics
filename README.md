# mutilomics_22q

This repository contains the code and data for the manuscript "".

# Introduction to 22q11.2 deletion syndrome




22q11.2 deletion syndrome (22q11.2DS) is the most common copy number variant (CNV) - associated syndrome, caused by the deletion of a small segment of chromosome 22.

> 22q11.2 deletion syndrome (22q11.2DS) is the most common chromosomal microdeletion disorder, estimated to result mainly from de novo non-homologous meiotic recombination events occurring in approximately 1 in every 1,000 fetuses. The first description in the English language of the constellation of findings now known to be due to this chromosomal difference was made in the 1960s in children with DiGeorge syndrome, who presented with the clinical triad of immunodeficiency, hypoparathyroidism and congenital heart disease. The syndrome is now known to have a heterogeneous presentation that includes multiple additional congenital anomalies and later-onset conditions, such as palatal, gastrointestinal and renal abnormalities, autoimmune disease, variable cognitive delays, behavioural phenotypes and psychiatric illness — all far extending the original description of DiGeorge syndrome. Management requires a multidisciplinary approach involving paediatrics, general medicine, surgery, psychiatry, psychology, interventional therapies (physical, occupational, speech, language and behavioural) and genetic counselling. Although common, lack of recognition of the condition and/or lack of familiarity with genetic testing methods, together with the wide variability of clinical presentation, delays diagnosis. Early diagnosis, preferably prenatally or neonatally, could improve outcomes, thus stressing the importance of universal screening. Equally important, 22q11.2DS has become a model for understanding rare and frequent congenital anomalies, medical conditions, psychiatric and developmental disorders, and may provide a platform to better understand these disorders while affording opportunities for translational strategies across the lifespan for both patients with 22q11.2DS and those with these associated features in the general population.
>
> 
![22q region gene](./images/22q_intro.png)


Here we list the genes in the 22q11.2 region:


![22q region gene](./images/22q_gene.png)
> Schematic representation of the 3-Mb 22q11.2 region that is commonly deleted in 22q11.2 deletion syndrome, including the four low copy repeats (LCR22s) that span this region (LCR22A, LCR22B, LCR22C and LCR22D). Common commercial probes for fluorescence in situ hybridization (FISH) are indicated (N25 and TUPLE). The protein-coding and selected non-coding (*) genes are indicated with respect to their relative position along chromosome 22 (Chr22). T-box 1 (TBX1; green box) is highlighted as the most widely studied gene within the 22q11.2 region. Mutations in this gene have resulted in conotruncal cardiac anomalies in animal models and humans. Known human disease-causing genes that map to the region are indicated in grey boxes. These include proline dehydrogenase 1 (PRODH; associated with type I hyperprolinaemia), solute carrier family 25 member 1 (SLC25A1; encoding the tricarboxylate transport protein and is associated with combined D-2- and L-2-hydroxyglutaric aciduria), platelet glycoprotein Ib β-polypeptide (GP1BB; associated with Bernard–Soulier syndrome), scavenger receptor class F member 2 (SCARF2; associated with Van den Ende–Gupta syndrome), synaptosomal-associated protein 29 kDa (SNAP29; associated with cerebral dysgenesis, neuropathy, ichthyosis and palmoplantar keratoderma (CEDNIK) syndrome), and leucine-zipper-like transcription regulator 1 (LZTR1; associated with schwannomatosis 2). Further details on the location of non-coding RNAs and pseudogenes in the 22q11.2 region may be found in Guna et al.89. Common 22q11.2 deletions are shown, with the typical 3-Mb deletion flanked by LCR22A and LCR22D (LCR22A– LCR22D) on top and the nested deletions, with their respective deletion sizes indicated below. Each of the deletions portrayed is flanked by a particular LCR22. Those rare deletions not mediated by LCRs are not shown. AIF3M, apoptosis-inducing factor mitochondrion-associated 3; ARVCF, armadillo repeat gene deleted in velocardiofacial syndrome; CDC45, cell division cycle 45; Cen, centromere; CLDN5, claudin 5; CLTCL1, clathrin heavy chain-like 1; COMT, catechol-O-methyltransferase; CRKL, v-crk avian sarcoma virus CT10 oncogene homologue-like; DGCR, DiGeorge syndrome critical region; GNB1L, guanine nucleotide-binding protein (G protein), β-polypeptide 1-like; GSC2, goosecoid homeobox 2; HIC2, hypermethylated in cancer 2; HIRA, histone cell cycle regulator; KLHL22, kelch-like family member 22; LINC00896, long intergenic non-protein-coding RNA 896; LOC101927859, serine/arginine repetitive matrix protein 2-like; CCDC188, coiled-coil domain-containing 188; LRRC74B, leucine-rich repeat-containing 74B; MED15, mediator complex subunit 15; mir, microRNA; MRPL40, mitochondrial ribosomal protein L40; P2RX6, purinergic receptor P2X ligand-gated ion channel 6; PI4KA, phosphatidylinositol 4-kinase catalytic-α; RANBP1, Ran-binding protein 1; RTN4R, reticulon 4 receptor; SEPT7, septin 7; SERPIND1, serpin peptidase inhibitor clade D (heparin co-factor) member 1; TANGO2, transport and golgi organization 2 homologue; THAP7, THAP domain-containing 7; TRMT2A, tRNA methyltransferase 2 homologue A; TSSK2, testis-specific serine kinase 2; TXNRD2, thioredoxin reductase 2; UFD1L, ubiquitin fusion degradation 1-like; USP41, ubiquitin-specific peptidase 41; ZDHHC8, zinc-finger DHHC-type-containing 8; ZNF74, zinc-finger protein 74.


**References:**

McDonald-McGinn, D., Sullivan, K., Marino, B. et al. 22q11.2 deletion syndrome. Nat Rev Dis Primers 1, 15071 (2015). https://doi.org/10.1038/nrdp.2015.71




# File structure
```bash
tree -L 2

.
├── README.md
├── analysis
│   ├── synaptosomes_bulkRNA_analysis
│   ├── synaptosomes_miRNA_anlaysis
│   ├── synaptosomes_miRNA_gene_network
│   └── synaptosomes_scRNA_analysis
├── data
│   ├── neruon_bulkRNA
│   ├── ref
│   ├── synaptosomes_bulkRNA
│   ├── synaptosomes_miRNA
│   └── synaptosomes_scRNA
├── images
│   ├── 22q_gene.png
│   └── 22q_intro.png
├── logs
└── source
    └── DEG_functions.R

15 directories, 4 files
```

