<h1><center>UM prognosis-related pathway analysis</h1></center>

![Flowchart](/Data/Flowchart.png)

# Data sets

Two bulk RNA-seq data for UM:
>* GSE176345: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176345.
>* TCGA data: https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Ocular%20melanomas%20(UVM)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443.

</br>

One scRNA-seq data for UM:
>* GSE139829: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139829.

# Procedures
1. The definition of survival information genes (SIGs).
> Corresponding code: `/Code/code1.R`

</br>

2. Single-cell data processing.
> Corresponding code: `/Code/code2.R`

</br>

3. Inferring the survival phenotype of cells.
> Corresponding code: `/Code/code3.R`

</br>

4. Pathway enrichment analysis.
> Corresponding code: `/Code/code4.R`

</br>

5. Identifying UM prognosis-related pathways.
>  Corresponding code: `/Code/code5.R`

</br>

6. Constructing a prognostic model for UM.
> Corresponding code: `/Code/code6.R`
The prognostic model of UM (*forest_model.Rdata*) and related data are stored in `./Data/`.

</br>

> The relevant plotting programs are stored in `./Code/Plot/`.
