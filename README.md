

A pipeline for creat Figure 2 from [Safra et al.](https://genome.cshlp.org/content/33/1/71.long)


##### Input files:

1. The read can downloaded from the NCBI Sequence Read Archive under BioProject accession ID: PRJNA788351
2. The Supplemental_Table_S1 from the Supplemental in [Safra et al.](https://genome.cshlp.org/content/33/1/71.long)
3. A table that map the SRR number to the namber in the Supplemental_Table
4. A germline files in fasta format.


##### Output files:

1. sample_germ-pass.fastq
2. sample_clone-pass.fastq
3. plot.pdf - Figure 2 from [Safra et al.](https://genome.cshlp.org/content/33/1/71.long)
4. Pipeairr_Figure2_Panelc.pdf .....


##### Pipeline container:

* Docker: immcantation/suite:4.3.0


##### Pipline processing steps:

1. fastqTofasta
2. Igblastn
3. MakeDb_igblast
4. CreateGermlines
5. change_names_germ_pass
6. DefineClones
7. change_names_clone_pass
8. crohn_analysis
9. crohn_plot_figure 2
10. pipeAIRR_figure2_panelC


#### Files:

* [Supplemental_Table_S1]()
* [map_table]()
