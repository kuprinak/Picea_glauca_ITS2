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
    B --> C@{shape: procs, label: "Trimmed reads (NCBI: PRJNA1335163)"};
    C --> Fa([FastQC]);
    Mu --> Mur@{shape: procs, label: "MultiQC_raw_reads.html, MultiQC_trimmed_reads.html"}; 
    C --> V([VSEARCH]);
    V --> T([UNITE v9.0 database]);
    T --> CV@{shape: procs, label: "Alaska_counts.txt"} --> R([R: Distribution_Growth.html, Metabarcoding.html, Metabarcoding_GAM_Growth-vs-RAF-diversity.html, Metabarcoding_Trophic-Guilds.html])
    T --> Ta@{shape: procs, label: "Alaska_taxonomy.txt"} --> R
    Ta --> F([FUNGuild v.1.1 database]);
    F --> Gu@{shape: procs, label: "Guilds and trophic modes of taxa"} --> R
    R --> alpha([Alpha diversity])
    R --> beta([Beta diversity])
    R --> Guild([Guild relative abundance])


Tree_data@{shape: procs, label: "Ring data, DBH"} --> R1([R: BAI_calculation.Rmd]) --> BAI@{shape: procs, label: "detrended BAI"} --> Sam_inf@{shape: procs, label: "Alaska_info.txt and Alaska_info_noR.txt"} --> R
Data@{shape: procs, label: "Tree age, height, soil pH"}   --> Sam_inf

```
