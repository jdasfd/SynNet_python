# SynNet: A python package of the SynNet pipeline

SynNet is a python package of the SynNet pipeline, which is rewritten from the original [SynNet-Pipeline](https://github.com/zhaotao1987/SynNet-Pipeline). We used the jcvi package (mcscan module) to build the micro-synteny networks.

## Installation

```bash
git clone https://github.com/jdasfd/SynNet_python.git
cd SynNet_python

conda env create -f environment.yml
conda activate synnet

pip install .
```

**Note**: If you want to develop/modify the code, use `pip install -e .` instead.

```bash
usage: synnet [-h] [-V] <command> ...

SynNet v0.1.0: Synteny Network Builder

Pipeline: gff2bed -> mcscan -> network -> cluster -> viz

positional arguments:
  <command>      Available commands
    gff2bed      Convert GFF3 to BED format
    mcscan       Chain-wise MCScan
    network      Build synteny network
    cluster      Cluster synteny network and filter
    viz          Visualize synteny network (interactive HTML)

options:
  -h, --help     show this help message and exit
  -V, --version  show program's version number and exit

Use 'synnet <command> --help' for command-specific help.
```

## Usage

### gff2bed

```bash
usage: synnet gff2bed [-h] -i INPUT [-o OUTPUT] [-s SPECIES_LIST] [--output-dir OUTPUT_DIR] [-t FEAT_TYPE] [-k ID_KEY]
                      [--min-length MIN_LENGTH] [-v]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input GFF3 file or directory [required]
  -o OUTPUT, --output OUTPUT
                        Output BED file (single file mode only)
  -s SPECIES_LIST, --species-list SPECIES_LIST
                        Species list file (one name per line)
  --output-dir OUTPUT_DIR
                        Output directory for BED files (default: same as input)
  -t FEAT_TYPE, --feat-type FEAT_TYPE
                        Feature type to extract (default: mRNA)
  -k ID_KEY, --id-key ID_KEY
                        Attribute key for gene ID (default: ID)
  --min-length MIN_LENGTH
                        Minimum feature length (default: 0)
  -v, --verbose         Verbose output
```

The gff2bed module is developed basically to convert gff3 to bed6 format. It is highly recommended to manipulate gff3 files on your own or prepare the bed6 file directly.

If your protein and annotation files could be accepted by jcvi, then everything will be fine.

```bash
# -i is the input dir: which will convert all gff3 to bed6 in the species.lst
synnet gff2bed -s species.lst -i ./
```

### mcscan

```bash
usage: synnet mcscan [-h] -s SPECIES_LIST [-i INPUT_DIR] [-o OUTPUT_DIR] [--cscore CSCORE] [--min-size MIN_SIZE]
                     [--cpus CPUS] [--no-intra] [--dry-run] [-v]

options:
  -h, --help            show this help message and exit
  -s SPECIES_LIST, --species-list SPECIES_LIST
                        Species list file [required]
  -i INPUT_DIR, --input-dir INPUT_DIR
                        Directory containing sequence and annotation files (default: current dir)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory for jcvi results (default: jcvi_output)
  --cscore CSCORE       C-score cutoff (default: 0.7)
  --min-size MIN_SIZE   Minimum anchors in a cluster (default: 4)
  --cpus CPUS           CPU cores (default: 4)
  --no-intra            Skip self synteny detection
  --dry-run             Print commands without executing
  -v, --verbose         Verbose output
```

The mcscan module is implemented by [`jcvi`](https://github.com/tanghaibao/jcvi). Better reading `jcvi` before you use this package.

The `species.lst` contains the list of the species. The inter-species mcscan will be carried out in the order of line-by-line in `species.lst`. This step will automatically create the `jcvi_output` directory.

`jcvi_output/`: intra- and inter-species anchors files

```bash
synnet mcscan -s species.lst -i ./ --cscore 0.99 --min-size 5
```

### network

```bash
usage: synnet network [-h] -s SPECIES_LIST [-i INPUT_DIR] [-o OUTPUT_DIR] [--no-lifted] [--min-score MIN_SCORE] [-v]

options:
  -h, --help            show this help message and exit
  -s SPECIES_LIST, --species-list SPECIES_LIST
                        Species list file [required]
  -i INPUT_DIR, --input-dir INPUT_DIR
                        Input directory containing anchors files (default: jcvi_output)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory for network files (default: network_output)
  --no-lifted           Exclude lifted alignments (rows with 'L' suffix in score)
  --min-score MIN_SCORE
                        Minimum score threshold (default: 0)
  -v, --verbose         Verbose output
```

This step will build the synteny network from `jcvi_output/*.lifted.anchors`. `.lifted.anchors` contains syntenic blocks. The network will be saved as the `network_output/Final_Network.tsv` file, containing the edges between nodes (with synteny relations).

```bash
synnet network -s species.lst
```

### cluster

```bash
usage: synnet cluster [-h] -i INPUT -s SPECIES_LIST -b BED_DIR [-o OUTPUT_DIR]
                      [--method {cc,louvain,infomap,label_prop,spectral}] [--k K] [--cluster-size CLUSTER_SIZE]
                      [--min-species MIN_SPECIES] [--gene-list GENE_LIST] [-v]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input network TSV file (from network step) [required]
  -s SPECIES_LIST, --species-list SPECIES_LIST
                        Species list file [required]
  -b BED_DIR, --bed-dir BED_DIR
                        Directory containing .bed files (for gene-species mapping) [required]
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory (default: network_output)
  --method {cc,louvain,infomap,label_prop,spectral}
                        Clustering method (default: cc)
  --k K                 Number of clusters for spectral method (default: 10)
  --cluster-size CLUSTER_SIZE
                        Minimum cluster size (default: 2)
  --min-species MIN_SPECIES
                        Minimum species count per cluster (default: 1)
  --gene-list GENE_LIST
                        Gene list file (same directory as species.lst) to filter clusters
  -v, --verbose         Verbose output
```

`network_output/` contains two major files:

`Filtered.clusters.tsv`: cluster ids and their related genes (format as the orthofinder orthogroups.tsv file).

`Filtered.cluster_summary.tsv`: cluster ids and their size, species count and detailed species composition.

This cluster module will also give out synteny clusters from the `network_output/Final_Network.tsv` file. The clusters will be saved as the `network_output/Clusters.synnet.tsv` file in default. `Clusters.synnet.tsv` will contain the cluster ID, the score, and two nodes with syntenic relations. If you want to extract some specific genes, use `--gene-list` to provide a gene list file (contains all genes you want to extract).

```bash
synnet cluster -i ./network_output/Final_Network.tsv -s species.lst -b ./
synnet cluster -i ./network_output/Final_Network.tsv -s species.lst -b ./ --gene-list ./gene.lst
```

### viz

```bash
usage: synnet viz [-h] -i INPUT -s SPECIES_LIST [-b BED_DIR] [-o OUTPUT_DIR] [-v]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Clusters.synnet.tsv file (from cluster command) [required]
  -s SPECIES_LIST, --species-list SPECIES_LIST
                        Species list file [required]
  -b BED_DIR, --bed-dir BED_DIR
                        Directory containing .bed files (default: current dir)
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory (default: same as input file)
  -v, --verbose         Verbose output
```

This viz module will provide a visualization result from `network_output/Clusters.synnet.tsv` file. The `.html` result could be seen and saved. But if there are too many clusters, the html may be slow or fail to load.

```bash
synnet viz -i network_output/Clusters.synnet.tsv -s species.lst
```
