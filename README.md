# TrH2_symbiosis

In this, you can read what I am doing during my research internship. 
The subject of my intership is to find out with which cell type of Trychoplax sp. H2 (TrH2) is associated with the symbionts. For this, the scRNA raw data of TrH2 were acquired from the NCBI. (Details in each code description.)The reference symbiont genomes were provided by Harald. (https://keeper.mpdl.mpg.de/d/f876ded599e94520ae27/)
The files provided were following:
- Symbiont reference genome .fasta, .gtf and .gff files
    - Grellia_incantans_TRH2_HVG.fasta
    - 6666666.192815.gbk.gff
        - This includes additional infos about the genome
        # This can be used to find out prevalent genes: replace gene ID of bc-feature matrix with the real gene name of the reference genome annotation
        # maybe try to summmarize/calculate the relative abundance with the host and symbiont gene counts
            - do it also for effective genes, 16S gene -> with statistics
            - find out the most and least abundant genes    
    - Grellia_incantans_TRH2_HVG.genbankCurated.gtf
        - contains transcript and gene ID of each transcript

    - RETAH2_CDS.fasta
    - Ruthmannia_H2.fasta
    - 

The scRNA-sequencing data is publicly available in the NCBI data base and can be found by accession numbers. Due to the big files, the data needs to be downloaded with sratoolkit (v3.1.1)

The cellranger was downloaded but I couldn't call the command inside the cellranger package, So Jinru ran it for me. (cellranger-8.0.1) He provided me following files from cellranger count runs:
- Command:
    cellranger count --id=SRR(accession no.) --create-bam=true --transcriptome=/gxfs_work/cau/sunzm503/20.HaraldLab/2.scTrichoplaxSymbiont/0.refGrellia/GrelliaIncantansGbk --fastqs=/gxfs_work/cau/sunzm503/20.HaraldLab/2.scTrichoplaxSymbiont/0.rawData/(path_to_sc_transcriptome_data)
-  For the scRNA samples with these (NCBI) accession numbers: 
    - SRR24887387
    - SRR24886407
    - SRR24886410
    - SRR24886411
    - SRR24886416
    - SRR24886417 
    - SRR24886420
    - SRR24886421
    - SRR24886422

- raw feature-barcode-matrix (gzipped): every barcide with at least one read contained, including background and cell-associated barcodes.
    - barcodes.tsv
    - features.tsv
    - matrix.mtx (row: features; column: barcodes)
- filtered feature-barcode-matrix (gzipped): Only cell-associated barcodes, background noises are filtered out
    - barcodes.tsv
    - features.tsv
    - matrix.mtx
- metrics_summary_csv
- web_summary
    - I summarized the web_sammary data of every SRR-samples in "TrH2_data_summary_csv" in "SEURAT" folder
- For details about inputs, running pipelines, ouputs and analysis, see: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis

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
        - SRR24886387 and SRR24886410, there are no estimated cell number and no reads mapped to the annotated reference genome despite of high read number
        - SRR24886416: this sample also has a low estimated cell number (2), but many reads mapped.   

        
    # If fails: how and why could this happen?
        reasons for no estimated cell and reads mapped despite of high read number might be:
        - abscence of intracellular symbiont in the sample
        - Existence of intracellular symbiont below the detection limit
        Reasons for low estimated cell number in general:
        - Different sample handling
        - Different condition of the examined Trichplax spp. individual, e.g. health, development, treatments etc. 
        -  Individual variance of intracellular symbiont associated with the Trichoplax individual
        - Of course, it is also possible that the host RNA is contained in the matrix due to the sequence similarity

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
        - Minimum barcode containing each features: 2 or 3

- Find out the actual genes and do statistics
    - Grellia_incantans_TRH2_HVG.genbankCurated.gtf.gff
        - This includes additional infos about the genome
        # This can be used to find out prevalent genes: replace gene ID of bc-feature matrix with the real gene name of the reference genome annotation
        # maybe try to summmarize/calculate the relative abundance with the host and symbiont gene counts
            - do it also for effective genes, 16S gene -> with statistics
            - find out the most and least abundant genes (only with positive cells)


SEURAT
SEURAT is a tool for single cell genomics analysis, particularly for the clustering analysis.  
For Seurat run, read SEURAT.ipynb (Saskianne/TrH2_symbiosis/SEURAT.ipynb) file. 




METACELL
Clustering and visualization 
For metacell run, read "Metacell.ipynb" file
For the visualization, read "metacell.r" file, mc2d download



