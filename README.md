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
    Mu --> Mur@{shape: procs, label: "MultiQC_report_trimmed_reads.html"}; 
    C --> V([VSEARCH]);
    V --> T([UNITE v9.0 database]);
    T --> CV@{shape: procs, label: "Alaska_counts.txt"} --> R([R: 2_Growth-pH.html, 3_Metabarcoding.html, 4_GLM: Growth-vs-RAF-diversity.html, 5_Trophic-Guilds.html])
    T --> Ta@{shape: procs, label: "Alaska_taxonomy.txt"} --> R
    Ta --> F([FUNGuild v.1.1 database]);
    F --> Gu@{shape: procs, label: "Guilds and trophic modes of taxa"} --> R
    R --> alpha([RAF Alpha diversity])
    R --> beta([RAF composition])
    R --> Guild([Guild relative abundance])


Tree_data@{shape: procs, label: "Ring data, DBH"} --> R1([R: 1_BAI_calculation.Rmd]) --> BAI@{shape: procs, label: "BAI, detBAI, rcsBAI"} --> Sam_inf@{shape: procs, label: "Alaska_info.txt and Alaska_info_noR.txt"} --> R
Data@{shape: procs, label: "Tree age, height, soil pH"}   --> Sam_inf

```
