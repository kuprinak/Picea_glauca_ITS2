
## Environment information:

Analysis was performed on the "High Performance Computing" (HPC) cluster of the University of Greifswald, 
which consists of 78 pure computing nodes @3.0 GHz with 64 GB RAM each and is managed by SLURM system with the installed:

conda v25.3.1 and Python v3.7.12

# 1. Nextera adapters removal in Trimmomatic

```{bash}
for f in /home/kuprinak/Diph/RawData/*_R1_001.fastq.gz

do
    n=${f%%_R1_001.fastq.gz} 
    trimmomatic PE -threads 96 $1 ${n}_R1_001.fastq.gz  ${n}_R2_001.fastq.gz \
    ${n}_R1_trimmed.fastq.gz ${n}_R1_unpaired.fastq.gz ${n}_R2_trimmed.fastq.gz \
    ${n}_R2_unpaired.fastq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:  \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60
done
```

# 2. Data Processing in VSEARCH


```{bash}
THREADS=24  # Select a desired number of threads for parallelization
REF=/home/kuprinak/Diph/DBs/ref.fasta # reference collection of fungi ITS  
CHIMREF=/home/kuprinak/Diph/DBs/uchime.fasta # reference collection of fungi ITS for chimera detection
```

### 2.1. Convert reference sequence databases into binary format

```{bash}
vsearch --makeudb_usearch $REF --dbmask none --output $REF.udb
vsearch --makeudb_usearch $CHIMREF --dbmask none --output $CHIMREF.udb
```

### 2.2. Process all raw FASTQ files one by one

>[!IMPORTANT]
>The commands below start a cycle for each pair of FASTQ files:  

* forward and reverse reads will be merged and saved in FASTQ format;
* low-quality reads will be removed and the remaining reads will be saved in FASTA format;
* reads will be dereplicated (only one copy of each unique sequence will be preserved, and information about its abundance in a sample will be added to its title).
  
```{bash}
for f in /home/kuprinak/Diph/DATA/RawTrimmomatic/*_R1_*.fastq.gz; do  # start a cycle to process every forward FASTQ file

    r=$(sed -e "s/_R1_/_R2_/" <<< $f)  # save a name of the reverse read as variable "r".

    s=$(awk -F'[_]' '{print $2}' <<< $f)  # save a name of the currently processed sample as "s". The program "awk" splits the file name by the "_" symbol and takes the second resulting element. For example, for file name "Run20211025Order279Sample003_2S1_L001_R1_001.fastq.gz" it is "2S1".

	echo
	echo ====================================
	echo Processing sample $s  # The command "Echo" prints into terminal the text that follows it.
	echo ====================================
```
### 2.3. Merging of forward and reverse reads

```{bash}
    vsearch --fastq_mergepairs $f \
    	--threads $THREADS \
        --reverse $r \
        --fastq_minovlen 150 \
        --fastq_maxdiffs 15 \
        --fastqout ./output/$s.merged.fastq \
        --fastq_eeout

    echo
    echo Calculate quality statistics
```
* fastq_mergepairs - the command itself, followed by the path to the file with forward reads (f)
* threads - desired number of parallel processes
* reverse - path to the file with reverse reads (r)
* fastq_minovlen - minimum allowed length of overlapping region of forward and reverse reads
* fastq_maxdiffs - maximum number of allowed mismatches in overlapping region
* fastqout - path for the output FASTQ file with merged reads
* fastq_eeout - add the estimated expected number of errors to the read name

### Calculation of the quality statistics for merged reads

```{bash}
    vsearch --fastq_eestats ./output/$s.merged.fastq \
        --output ./output/$s.stats
    echo
    echo Quality filtering
```

### 2.4. Remove low-quality reads from the samples
    
```{bash}
    vsearch --fastq_filter ./output/$s.merged.fastq \
        --fastq_maxee 0.5 \
        --fastq_minlen 250 \
        --fastq_maxns 0 \
        --fastaout ./output/$s.filtered.fasta \
        --fasta_width 0    

    echo
    echo Dereplicate at sample level and relabel with sample_n
```
  * fastq_filter - the command itself, followed by the path to the file with merged sequences.
  * fastq_maxee - maximum allowed expected number of errors in a sequence
  * fastq_minlen - minimum allowed sequence length
  * fastq_maxns - max allowed number of uncertain bases (Ns)
  * fastaout - path for the output FASTA file
  * fasta_width - length of the lines for the output FASTA file. 0 = no line splitting
  * fastq_stripleft/stripright - how many characters to strip from both ends to remove primer regions.
    
 ### 2.5. Sequence dereplication
 
```{bash}
    vsearch --derep_fulllength ./output/$s.filtered.fasta \
        --strand both \
        --output ./output/$s.derep.fasta \
        --sizeout \
        --uc ./output/$s.derep.uc \
        --relabel $s. \
        --fasta_width 0

done
```
  * derep_fulllength - the command itself, followed by the path to the file with quality-filtered sequences
  * strand - "both" to merge complementary sequences
  * output - path for the output FASTA file
  * sizeout - preserve the original abundance of sequence as part of the sequence name
  * uc - path to the output "map" of sequences in UC format
  * relabel - give a new name (here - sample name) to each sequence followed by 1, 2, 3 etc.
    
>[!IMPORTANT]
> END of the cycle
    
### 2.6. Pool all samples together and dereplicate globally
 
#### Merge all samples

```{bash}
cat ./output/*.derep.fasta > ./output/all.fasta  # Pool sequences from all samples together
```
#### Dereplicate across samples and remove singletons

* minuniquesize - minimum allowed abundance of a sequence. If 2, all singletons will be removed. 
* sizein - sequence abundance information stored in the sequence name will be taken into account
  
```{bash}
vsearch --derep_fulllength ./output/all.fasta \
    --strand both \
    --minuniquesize 10 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc ./output/all.derep.uc \
    --output ./output/all.derep.fasta

echo Unique non-singleton sequences: $(grep -c "^>" ./output/all.derep.fasta)
```
* minuniquesize - minimum allowed abundance of sequence. If 2, all singletons will be removed. 
* sizein - sequence abundance information stored in the sequence name will be taken into account
  
### 2.7. Preclustering

The resulting sequences are sorted by abundance (most abundant go first).

```{bash}
vsearch --cluster_size ./output/all.derep.fasta \
    --id 0.99 \
    --threads $THREADS \
    --strand both \
    --fasta_width 0 \
    --sizein \
    --sizeout \
    --uc ./output/preclustered.uc \
    --centroids ./output/preclustered.fasta

echo Unique sequences after preclustering: $(grep -c "^>" ./output/preclustered.fasta)
```

### 2.8. Chimera detection

Two chimera detection algorithms will be applied:
* De novo chimera detection - based on comparison of analyzed sequences with each other;
* Reference-based detection - based on comparison of analyzed sequences with reference sequences.

```{bash}
vsearch --uchime3_denovo ./output/preclustered.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --qmask none \
    --nonchimeras ./output/denovo.nonchimeras.fasta \

echo sequences after de novo chimera detection: $(grep -c "^>" ./output/denovo.nonchimeras.fasta)

vsearch --uchime_ref ./output/denovo.nonchimeras.fasta \
	--threads $THREADS \
    --db $CHIMREF \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --qmask none \
    --dbmask none \
    --nonchimeras ./output/nonchimeras.fasta

echo Unique sequences after reference-based chimera detection: $(grep -c "^>" ./output/nonchimeras.fasta)

```

### 2.9. Clustering sequences to OTUs (greedy clustering)

```{bash}
vsearch --cluster_size ./output/nonchimeras.fasta \
    --relabel OTU \
    --id 0.97 \
    --threads $THREADS \
    --strand both \
    --fasta_width 0 \
    --sizein \
    --sizeout \
    --uc ./output/OTUs.uc \
    --centroids ./output/OTUs.fasta

echo resulting OTUs: $(grep -c "^>" ./output/OTUs.fasta)

```

### 2.10. Map sequences to OTUs by searching

```{bash}
vsearch --usearch_global ./output/all.fasta \
    --threads $THREADS \
    --db ./output/OTUs.fasta \
    --id 0.97 \
    --strand both \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --qmask none \
    --dbmask none \
    --otutabout ./output/OTUtable.txt

# echo Sort OTU table numerically (by decreasing OTU number)
sort -k1.4n ./output/OTUtable.txt > ./output/OTUtable.sorted.txt # -k1.4 means sorting by the 4th symbol of the first column
```
* usearch_global - the command itself, followed by a path to quality-filtered sequences pooled from all samples (but not yet dereplicated globally)
* db - path to the file with OTUs
* id - sequence similarity threshold (from 0 to 1)
* otutabout - path for the output OTU table in TSV format
  
### 2.11. Searching OTUs across reference sequence database

 OTUs will be compared with reference sequences and results written into a table:   | OTU | Reference | Similarity |  

#### Find the best matches of OTU centroids with the fungi database (UNITE DB)
##### Search OTU centroids across UNITE database

```{bash}
vsearch --usearch_global ./output/OTUs.fasta \
    --threads $THREADS \
    --db $REF.udb \
    --qmask none \
    --dbmask none \
    --id 0.95 \
    --top_hits_only \
    --mincols 250 \
    --maxaccepts 100 \
    --maxrejects 100 \
    --userout ./output/OTUs.fungi.tsv \
    --userfields query+target+id \
    --strand both \
    --alnout ./output/OTUs.fungi.alignment.fasta
```

# 3. FUNGuld search

```{bash}
python3 Guilds_v1.1.py -otu FunGuild_OTU.txt -db fungi 
```
Output: "..."
