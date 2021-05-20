# Genomic Data Science Project: RNA-seq for Studying Brain Development 
  
In this [project](https://www.coursera.org/learn/genomic-data-science-project) I will replicate the data analysis of the RNA-seq experiment conducted in ([Jaffe et al., 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4281298/)).
   
The aim of this experiment is to study the development of the human brain across the lifespan. The data generated measured gene expression from different individuals across the human lifespan. The authors studied 6 different age groups: fetal (<0 years), infant (0-1 years), child (1-10 years), adolescent (10-20 years), adult (20-50 years) and old, (50+ years), and each age group had 6 indiividuals. They were looking for genes that showed patterns of expression that changed between the 6 age groups.
  
Here, seeking the manageability of the data, I only consider 2 out of 6 age groups: fetal (<0 years) and adult (20-50 years) age groups, and each age group had 3 individuals only.
  
PDF reports describing the, step by step, analysis are provided inl:
  * [Report-P1-Getting_the_Data.pdf](Report-P1-Getting_the_Data.pdf)
  * [Report-P2-Aligning_the_Reads.pdf](Report-P2-Aligning_the_Reads.pdf)
  * [Report-P3-Counting_the_Reads.pdf](Report-P3-Counting_the_Reads.pdf)
  * [Report-P4-Exploratory_Analysis.pdf](Report-P4-Exploratory_Analysis.pdf)
  * [Report-P5-Differential_Expression_Analysis.pdf](Report-P5-Differential_Expression_Analysis.pdf)
  * [Report-P6-Epigenetics_Analysis.pdf](Report-P6-Epigenetics_Analysis.pdf)
 
The corresponding Markdown files are provided in [Notebooks-RMarkdown](Notebooks-RMarkdown).  

Moreover,the bash and R scripts used in the analysis are provided in [bash-scripts](bash-scripts) and [R-scripts](R-scripts) respectively.
  
### Reference:   

[1] Jaffe AE, Shin J, Collado-Torres L, et al. [Developmental regulation of human cortex transcription and its clinical relevance at base resolution](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4281298/). Nature neuroscience. 2015;18(1):154-161.
