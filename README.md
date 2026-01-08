# snakemake pipeline for sciT analysis

This repository contains the scripts and instructions you will need to perform a sciT analysis.

This snakemake workflow represents an adaptation of an existing in-house snakemake workflow, which allows
it to preprocess single-cell transcriptomics data that was generated using a single-cell combinatorial
indexing approach (sciT).

## initial setup

### Install miniforge3

The first time you run the sciT-snakemake pipeline, you have to prepare an environment that can prepare and run the neccesary code.

You can use conda for this. I recommend you to use *Miniforge3* as a conda package manager. It can be installed like this:

```bash
# Use wget to download the installer
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
# Or use curl to download the installer
#curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

# Run the installer
bash Miniforge3-$(uname)-$(uname -m).sh
```

> [!NOTE]
> Using the installer like above should run it in interactive mode, allowing you to indicate the install directory for Miniforge3.

> [!NOTE]
> The installer will prompt you to initialize conda with your shell. Select yes if you want conda to activate automatically.

> [!IMPORTANT]
> If installing miniforge3 in an HPC, do **NOT** install it in your home directory. It will quickly run out of storage space for your environments.
>

### Create base conda environment

Now you can create the base conda environment that will be used for running the sciT workflow.

After conda is activated, run the following code to create the base environment.

```bash
conda create -c conda-forge -c bioconda -n sciT_snakemake -y snakemake samtools p7zip gcc zlib openjdk snakemake-executor-plugin-slurm mamba wget git
```

### Check reference repositories

For the sciT workflow to work, you need proper reference files. These can be placed in a `reference_repositories` directory (When executing the sciT workflow we will link to the exact location of the `reference_repositories` directory, so it can be anywhere you like). This should include:

- **contaminants:** A directory with *bowtie2* genome indexes and fasta files for potential contaminants to test for using *fastq-screen*. for example:

  ```txt
  |____reference_repositiories
  | |____contaminants
  | | |____cutibacterium.1.bt2
  | | |____cutibacterium.2.bt2
  | | |____cutibacterium.3.bt2
  | | |____cutibacterium.4.bt2
  | | |____cutibacterium.rev.2.bt2
  | | |____cutibacterium.rev.1.bt2
  | | |____cutibacterium.fasta
  ```

- Other directories containing genome files and gtf annotations for species you want to map your read onto. For example, here is a directory for **mm39 references**:

  ```txt
  |____reference_repositiories
  | |____mm39
  | | |____annotations.gtf
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.1.bt2
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.2.bt2
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.3.bt2
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.4.bt2
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.rev.1.bt2
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.rev.2.bt2
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.fa
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.fai
  | | |____Mus_musculus.GRCm39.dna_sm.primary_assembly.dict
  | | |____mgp_REL2021_snps.vcf.gz.csi
  | | |____mgp_REL2021_snps.vcf.gz.tbi
  | | |____mgp_REL2021_snps.vcf.gz
  ```

  Here, there are

  - **bowtie index files:** for including *mm39* in the fastq-screen report.
  - **fasta file:** for building the index for the [STAR aligner](https://pubmed.ncbi.nlm.nih.gov/23104886/).
  - **annotations.gtf:** For gene counting by [HTseq](https://htseq.readthedocs.io/en/latest/).

> [!NOTE]
> For the references to be included in the *fastq-screen* report, they need to be included in the `fastq_screen.conf` file.
> 
>     DATABASE	Mouse	reference_repositories/mm39/Mus_musculus.GRCm39.dna_sm.primary_assembly
>     DATABASE	Cutibacterium_Acnes	reference_repositories/contaminants/cutibacterium
  
When you run the sciT worklfow, you will select *one reference fasta file* and one *GTF annotation file*.
This will be selected in the *config.yaml* file (more on that later).

### Install sciT workflow packages

In order to make the sciT workflow work, we need to install some additional packages. These are located in the `software_repositories` directory in this github repository. To install these packages in the **sciT_snakemake** conda environment, first activate the environment:

```bash
conda activate sciT_snakemake
```

Clone this github repository

```bash
# clone repository
git clone https://github.com/Suirotras/SciT_snakemake.git
# move into repository
cd SciT_snakemake
```

Then install the packages into the conda environment using **pip**.

```bash
pip install ./software_repositories/* --upgrade
```

This installation only needs to be performed once, after which this **sciT_snakemake** environment can be used for different sciT workflows.

You can now delete the cloned github repository, as you only used it here to install the sciT workflow packages.

## Execution of sciT workflow

### Activate the base environment.

```bash
conda activate sciT_snakemake
```

### Clone repository

Clone this github repository to the location where you want to preprocess the sciT fastq.gz files.

```bash
cd path/to/desired/location
# clone repository
git clone https://github.com/Suirotras/SciT_snakemake.git
# move into repository
cd SciT_snakemake
cd sciT
```

> [!IMPORTANT]
> If this is the FIRST time you are running this workflow, you NEED to first install the packages in the `software_repositories` directory.
> Check out the [**Install sciT workflow packages**](#install-scit-workflow-packages) section first before continuing.

### Create a symlink to the reference repository directory

You now need to tell where the `reference_repositories` live, and subsequently create a symbolic link to that directory. If you do not do this, then the sciT workflow will not be able to find the references.

```bash
ln -s ../reference_repositories reference_repositories
```

### Modify config.yaml

You will need to configure the configuration file that tells the sciT workflow which options to use and which files to process. Here is an example of a correct `config.yaml`:

```yaml
reference:
  # Human Hg38:
  # fasta: reference_repositories/HG38.p14/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
  # annotation_gtf: reference_repositories/HG38.p14/annotations.gtf
  # For mouse mm10 use:
  fasta: reference_repositories/mm10_pool_lab/Mus_musculus.GRCm38.dna.primary_assembly.with_ERCC.fa
  annotation_gtf: reference_repositories/mm10_pool_lab/annotations.gtf
  # hybrid: # Comment or remove this section to not perform allele specific analysis
  #   multisample_vcf_file: reference_repositories/mm39/mgp_REL2021_snps.vcf.gz
  #   alleleA: 129S1_SvImJ # Corresponds to a sample name in the VCF file above
  #   alleleB: CAST_EiJ

QC: # The quality control thresholds used for the quality control report
  transcriptome_count_threshold: 0
  transcriptome_genes_threshold: 1000
  report_top_samples: 100 # Number of best cells/samples to include in the report, based on the number of unique sciT reads
  report_bottom_samples: 100 # Number of worst cells to include in the report, based on the number of unique sciT reads

sciT:
  cell_min_transcriptome_count: 25 # Minimum number of transcriptome molecules for a cell to be counted

conditions:
  attributes: # Attributes used to determine the condition (label) of each sample
  - number
  - library
  colors: # Predetermined colors for attribute values:
    number:
      1K: 0000ff
      10K: abc8a4
    library:
      i31: ffcc00

# Add library data here by creating a libraries section, see documentation for more details
libraries:
  i31-sample1-1K:
    type: sciT
    demultiplexer: sciT
    transcriptome-paired-end-fastq-files:
      R1:
        - i31-sample1_L001_R1_001.fastq.gz
        - i31-sample1_L002_R1_001.fastq.gz
      R2:
        - i31-sample1_L001_R2_001.fastq.gz
        - i31-sample1_L002_R2_001.fastq.gz
    library-meta:
      number: 1K
      library: i31
    cell-meta:
      primary-index: sciT
      meta: {}
```

`config.yaml` consists of the following fields:

- **reference**

  - **fasta:** genome file from which STAR index will be created. Reads will be mapped to this STAR index.
  - **annotation_gtf:** GTF file that will be used to count transcriptomic features.
  - **hybrid:** Uncomment to perform an allele-specific analysis
    - **multisample_vcf_file:** vcf file needed for allele-specific analysis.
                                Uncomment to perform an allele-specific analysis.
    - **alleleA:** sample name in vcf file of allele A. Uncomment to perform an allele-specific analysis.
    - **alleleB:** sample name in vcf file of allele B. Uncomment to perform an allele-specific analysis.

- **QC:** Settings that will influence the QC report

- **sciT**
  - **cell_min_transcriptome_count:** A cell filter. Will filter out cells with less counts than the number indicated here.

- **conditions**
  - **attributes:** These attributes, in conjunction with the values in the *libraries/library-meta* fields, are used to select colors that will be used in
                    the quality report. For example, a **seaborn clustermap** that is generated will use these colors.

- **libraries:** Here you put the information about individual libraries that will be processed.

  - **i31-sample1-1K:** Should be changed to the library identifier of your choosing.
    - **type:** NEEDS to be **sciT** for every library entry.
    - **demultiplexer:** NEEDS to be **sciT** for every library entry.
    - **transcriptome-paired-end-fastq-files:** Indicate the paired end fastq-files that represent this library and should be processed.
    - **library-meta:** This indicates the library-wide metadata that defines this library. This metadata will be added to the generated
                        loom files that the sciT workflow will output.

### Modify fastq_screen.conf

In order for **fastq-screen** to work correctly, make sure the repositories are correctly listed in the `fastq_screen.conf`. Check this for both the
contaminants and references you are using for read mapping. This allows **fastq-screen** to find the correct bowtie genome indexes.

### Execute the sciT workflow

With all this set up, you can start the sciT workflow in order to process the sciT data.

#### for usage locally (i.e. your own laptop)

Adjust the core count and memory size

```bash
snakemake --cores 32 --resources mem_mb=62000 -k -p --rerun-incomplete --nt --use-conda --conda-prefix environments
```

#### for usage on SLURM cluster

The `--cores 60` and `mem_mb=100000` might not be the most ideal settings, You can probably change it to lower values.

```bash
snakemake --executor slurm --cores 60 --resources mem_mb=100000 -k -p --rerun-incomplete --jobs 20 --restart-times 3 --use-conda --conda-prefix environments
```

> [!NOTE]
> Use the `screen` utility to run the sciT workflow. This makes sure the sciT workflow can keep running when you to disconnect from the HPC.


## Issues

If there are issues, the documentation of the original in-house workflow might be helpful resolving the issue.