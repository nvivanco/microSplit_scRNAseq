# Notes

## Repo Files and Directories

- **`scripts/`**: containing bash and SLURM scripts for scRNA analysis, database download, and other utilities  
  - **`download_refgenome.sh`**: script for downloading reference genomes, transcriptomes, and GTFs  
  - **`download_sra.sbatch`**: SLURM script to download SRA samples from NCBI  
  - **`fastqc.sbatch`**: SLURM script to run FastQC  
  - **`star_index.sh`**: bash script to index reference genome with STAR  
  - **`star_align.sbatch`**: SLURM script to map reads to the indexed reference genome and run STARSolo  

Other files in the repository:  
- **`setup.sh`**: bash script to set up container and default paths  
- **`Dockerfile`**: Docker image to build container  
- **`README.md`**: project README  
- **`Notes.md`**: running documentation  

---

## Usage

- **`download_refgenome.sh`**  
  ```bash
  ./download_refgenome.sh <output_directory_path>
  ```

- **`download_sra.sbatch`**  
  ```bash
  sbatch download_sra.sh SRP SRP266243 $RAW_DATA/sra
  ```

- **`fastqc.sbatch`**  
  ```bash
  sbatch fastqc.sbatch $RAW_DATA/sra/SRP266243/SRR11940660 \
                       $PROCESSED_DATA/SRP266243/SRR11940660/fastqc_outs
  ```

- **`star_index.sh`**  
  ```bash
  ./star_index.sh <input_fasta> <input_gtf> <output_folder>
  # Example:
  ./star_index.sh $ECOLI_REF_GENOME $ECOLI_GTF \
                  $PROCESSED_DATA/ref_index/STAR_ecoli_genome_index
  ```

- **`star_align.sbatch`**  
  ```bash
  # Usage:
  sbatch star_align.sbatch \
         $PROCESSED_DATA/ref_index/STAR_bsub_genome_index \
         $RAW_DATA/sra/SRP266243/SRR11940660/SRR11940660_1.fastq \
         $RAW_DATA/sra/SRP266243/SRR11940660/SRR11940660_2.fastq \
         $RAW_DATA/barcodes/SRP266243/barcodes_r2_r3_solo.txt \
         $RAW_DATA/barcodes/SRP266243/barcodes_r1_solo.txt \
         $PROCESSED_DATA/SRP266243/SRR11940660/STAR_Bsub_ref_outs
  ```

---

## Data Organization

**`scRNAseq/`**  
- **`envs/`**: directory containing environment files  
  - `scrna_latest.sif`: singularity image file  

- **`processed_data/`** (`$PROCESSED_DATA`): contains analysis outputs  
  - **`ref_index/`**: indexed reference genomes  
    - `STAR_alt_bsub_genome_index`: `star_index.sh` output for reference genome of B.subtilis PY79  
    - `STAR_bsub_genome_index`: `star_index.sh` output for reference genome of B.subtilis strain 168  
    - `STAR_ecoli_genome_index`: `star_index.sh` output for reference genome of E.coli K-12 substr. MG1655  
  - **`SRP266243/`**: outputs for SRP266243  
    - **`SRR11940660/`**: outputs for SRR11940660 runs  
      - `fastqc_outs`: `fastqc.sbatch` output  
      - `STAR_alt_bsub_ref_outs`: STAR outputs (aligning reads to `STAR_alt_bsub_genome_index`)  
      - `STAR_Bsub_ref_outs`: STAR outputs (aligning reads to `STAR_bsub_genome_index`)  

- **`raw_data/`** (`$RAW_DATA`): containing raw data files  
  - **`barcodes/`**: barcodes from [Supplementary information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41596-024-01007-w/MediaObjects/41596_2024_1007_MOESM2_ESM.xlsx) from this [paper](https://www.nature.com/articles/nmeth.4220)  
    - `SRP266243/`  
  - **`references/`**  
    - `bsubtilis`: reference genomes/transcriptomes/gtf files of B.subtilis strain 168 and PY79  
    - `ecoli`: reference genomes/transcriptomes/gtf files of E.coli K-12 substr. MG1655  
  - **`sra/`**  
    - `SRP266243/`: this project has 3 runs  
      - `SRR11940660`: B.subtilis rep 1  
      - `SRR11940661`: B.subtilis rep 2  
      - `SRR11940662`: B.subtilis + E.coli  

---

## Methods

0. Run `setup.sh` script to set up environment variables and singularity run  

1. Prepare necessary files in `$RAW_DATA`  

Download the reference files:  
```bash
./download_refgenome.sh $RAW_DATA
```

Download SRA reads:  
```bash
sbatch download_sra.sh SRP SRP266243 $RAW_DATA/sra
```

Copy barcodes from Supplementary information to text files in `$RAW_DATA/barcodes`  

2. Convert file formats  

TODO: add gff -> gtf conversion here  

3. Run `fastqc.sbatch` on SRR11940660 (B. subtilis experiment):  
```bash
sbatch fastqc.sbatch $RAW_DATA/sra/SRP266243/SRR11940660 \
                     $PROCESSED_DATA/SRP266243/SRR11940660/fastqc_outs
```

4. Index reference genomes  

E.coli:  
```bash
./star_index.sh $ECOLI_REF_GENOME $ECOLI_GTF \
                $PROCESSED_DATA/ref_index/STAR_ecoli_genome_index
```

B.subtilis strain 168:  
```bash
./star_index.sh $BSUB_REF_GENOME $BSUB_GTF \
                $PROCESSED_DATA/ref_index/STAR_bsub_genome_index
```

B.subtilis strain PY79:  
```bash
./star_index.sh $BSUB_PY79_REF_GENOME $BSUB_PY79_GTF \
                $PROCESSED_DATA/ref_index/STAR_alt_bsub_genome_index
```

5. STARsolo align  

B.subtilis strain 168:  
```bash
sbatch star_align.sbatch \
       $PROCESSED_DATA/ref_index/STAR_bsub_genome_index \
       $RAW_DATA/sra/SRP266243/SRR11940660/SRR11940660_1.fastq \
       $RAW_DATA/sra/SRP266243/SRR11940660/SRR11940660_2.fastq \
       $RAW_DATA/barcodes/SRP266243/barcodes_r2_r3_solo.txt \
       $RAW_DATA/barcodes/SRP266243/barcodes_r1_solo.txt \
       $PROCESSED_DATA/SRP266243/SRR11940660/STAR_Bsub_ref_outs
```

B.subtilis strain PY79:  
```bash
sbatch star_align.sbatch \
       $PROCESSED_DATA/ref_index/STAR_alt_bsub_genome_index \
       $RAW_DATA/sra/SRP266243/SRR11940660/SRR11940660_1.fastq \
       $RAW_DATA/sra/SRP266243/SRR11940660/SRR11940660_2.fastq \
       $RAW_DATA/barcodes/SRP266243/barcodes_r2_r3_solo.txt \
       $RAW_DATA/barcodes/SRP266243/barcodes_r1_solo.txt \
       $PROCESSED_DATA/SRP266243/SRR11940660/STAR_alt_bsub_ref_outs
```