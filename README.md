# _Picea_ _glauca_ metabarcoding vs growth
A study examining the connection between root-associated fungi and the growth performance of _Picea glauca_.

> [!TIP]
> To view any HTML file, please download it


## Workflow:

```mermaid
flowchart TB
    A@{shape: procs, label: "Illumina raw reads"} --> B([Trimmomatic]);
    A --> Fa;
    Fa --> Mu([MultiQC]);
    B --> C@{shape: procs, label: "Trimmed reads (NCBI: PRJNA1185013)"};
    C --> Fa([FastQC]);
    Mu --> Mur@{shape: procs, label: "MultiQC_raw_reads.html, MultiQC_trimmed_reads.html"}; 
    C --> V([VSEARCH]);
    V --> T([UNITE v9.0 database]);
    T --> CV@{shape: procs, label: "Raw_Counts.txt"};
    T --> Ta@{shape: procs, label: "Taxonomy.txt"};
    Ta --> F([FUNGuild v.1.1 database]);
    F --> Gu@{shape: procs, label: "Guilds and trophic modes of taxa"};
    Gu --> rs([rstatix])
    CV --> D([DESeq2]);
    D --> N@{shape: procs, label: "Normalized and filtered read counts"};
    Ta --> Mi;
    Ta --> Me
    Ta --> Ph
    N --> Ph([phyloseq])
     N --> Mi([microeco])
    N --> Me([metacoder])
    Ta --> iN([iNEXT]);

```
