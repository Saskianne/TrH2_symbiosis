# TrH2_symbiosis

In this, you can read what I am doing during my research internship. 
The subject of my intership is to find out with which cell type of Trychoplax sp. H2 (TrH2) is associated with the symbionts. For this, the scRNA raw data of TrH2 were acquired from the NCBI. (Details in each code description.)The reference symbiont genome were provided by Harald. 
The files provided were following:
- Symbiont reference genome .fasta and .gtf files
    - Grellia_incantans_TRH2_HVG.fasta
    - Grellia_incantans_TRH2_HVG.genbankCurated.gtf.gff
        - This includes additional infos about the genome
        # This can be used to find out prevalent genes: replace gene ID of bc-feature matrix with the real gene name of the reference genome annotation
        # maybe try to summmarize/calculate the relative abundance with the host and symbiont gene counts
            - do it also for effective genes, 16S gene -> with statistics
            - find out the most and least abundant genes
    - RETAH2_CDS.fasta
    - Ruthmannia_H2.fasta
    - 
- 

The cellranger was downloaded but I couldn't call the command inside the cellranger package, So Jinru ran it for me. (cellranger-8.0.1) He provided me following files from cellranger count runs:
- Command:
    cellranger count --id=SRR24886387 --create-bam=true --transcriptome=/gxfs_work/cau/sunzm503/20.HaraldLab/2.scTrichoplaxSymbiont/0.refGrellia/GrelliaIncantansGbk --fastqs=/gxfs_work/cau/sunzm503/20.HaraldLab/2.scTrichoplaxSymbiont/0.rawData/SRR24886387
-  For the scRNA samples with these SRR numbers: 
    - SRR24887387
    - SRR24886407
    - SRR24886410
    - SRR24886411
    - SRR24886416
    - SRR24886417 
    - SRR24886420
    - SRR24886421
    - SRR24886422

- raw feature-barcode-matrix (gzipped)
    - barcodes.tsv
    - features.tsv
    - matrix.mtx
- filtered feature-barcode-matrix (gzipped)
    - barcodes.tsv
    - features.tsv
    - matrix.mtx
- metrics_summary_csv
- web_summary
    - I summarized the web_sammary data of every SRR-samples in "TrH2_data_summary_csv" in "SEURAT" folder

The matrices must be filtered by several criteria:
- scRNA sequencing quality parameters:
    - estimated number of cells
    - Reads
        - total read no. 
        - mean reads per cell
    - Genes
        - total genes detected
        - median genes per cell
    - Median UMI counts per cell
    # According to these criteria, certain samples were excluded from further analysis
    # If fails: how and why could this happen?

- Endosymbiont-relevant criteria:
    - genes mapped total to endosymbiont refererence
        - in total
        - median no. genes per cell
    - confidence of the reads mapped
    # narrow down the feature_bc_matrix
    - Pass/fail to exceed certain threshold
        - minimum feature count per barcode (gene per cell), maybe 3?
        # must be reasoned why
        - only barcodes with at least 1 features

- Find out the actual genes and do statistics
    - Grellia_incantans_TRH2_HVG.genbankCurated.gtf.gff
        - This includes additional infos about the genome
        # This can be used to find out prevalent genes: replace gene ID of bc-feature matrix with the real gene name of the reference genome annotation
        # maybe try to summmarize/calculate the relative abundance with the host and symbiont gene counts
            - do it also for effective genes, 16S gene -> with statistics
            - find out the most and least abundant genes (only with positive cells)


SEURAT

The SEURAT package was also resistant against running. Still I got it.
See "seurat_installation_script.r"
For Seurat run, read Seurat.ipynb file


METACELL
Clustering and visualization 
For metacell run, read "Metacell.ipynb" file
For the visualization, read "metacell.r" file, mc2d download



