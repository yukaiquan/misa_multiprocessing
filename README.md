# Genome-wide multi-process chromosome-level **SSR** search.

## 1 Installation

 *linux-64*    yukaiquan 1962568272@qq.com

Required software tools (third party, please install them according to instructions provided at respective links; e.g. by a conda environment):

```shell
conda install -c bioconda samtools
conda install -c conda-forge pandas
conda install -c bioconda pybedtools
```

## 2 Usage

### 2.1 help:

```shell
python misa_multiprocessing.py -h
usage: python misa_multiprocessing.py [-h] [--version] -g GENOME [-c CHR] [-t THREADS]

whole genome ssr find and primer design. Author:yukaiquan; Email:1962568272@qq.com

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  -g GENOME, --genome GENOME
                        Fatsa genome, should include all sequences from genome file.
  -c CHR, --chr CHR     Chromosome name, should be in the format of chr1A,chr2A,chr3A...(default is all chromosomes).
  -t THREADS, --threads THREADS
                        The number of threads (default is 1).

Example:
python kq_misa.py -g SFS.fasta -t 1\
```

### 2.2 example：

#### work folder

```shell
ls
Avena_eriantha.fasta.gz  Avena_eriantha.fasta.gz.fai  Avena_eriantha.fasta.gz.gzi  misa.ini  misa_multiprocessing.py  misa.pl
```

I have a genome file with 7 chromosomes.

```shell
cat Avena_eriantha.fasta.gz.fai
AE1     588203704       5       60      61
AE2     555382095       598007110       60      61
AE3     551069265       1162645579      60      61
AE4     534821622       1722899337      60      61
AE5     529955746       2266634658      60      61
AE6     459891171       2805423005      60      61
AE7     455803086       3272979034      60      61
chrUn   102925192       3736378846      60      61
```

I can 7 processes in parallel

```shell
python misa_multiprocessing.py -g Avena_atlantica.fasta.gz -t 7
```

**-g genomefile（gz/fasta）**

**-t process (default=1)**

**-c chr name （default=all）**

#### result

```shell
ls
Avena_atlantica.fasta.gz_second.tsv        chr2A.fasta_second.fasta.statistics  chr5A.fasta_second.fasta.statistics
Avena_atlantica.fasta.gz_second_all.fasta  chr3A.fasta_second.fasta.statistics  chr6A.fasta_second.fasta.statistics
chr1A.fasta_second.fasta.statistics        chr4A.fasta_second.fasta.statistics  chr7A.fasta_second.fasta.statistics
##########################################################################################
head -n 5 Avena_atlantica.fasta.gz_second.tsv
ID      SSR nr. SSR type        SSR     size    start   end
chr1A:42775-43110       1       c       (AC)6(AT)12     36      150     185
chr1A:174688-175102     1       c       (TTG)5ttctgattctgattttgcaggaactgtgtcatcatctgcattagttgagcttgtcctgtcaacatctgctgtatctcgggagtga(TGT)5       115     150     264
chr1A:182555-182866     1       p2      (GA)6   12      150     161
chr1A:189505-189816     1       p2      (GA)6   12      150     161
##########################################################################################
head -n 5 Avena_atlantica.fasta.gz_second_all.fasta
>chr1A:42775-43110
AAATTTAATATTTTAGAATAATTTACATTGCTTACTCAAGCAAAGGGACTCAATCTCGTGAATGTTTAAACGCTACATGCTCAAGCTGTTGTCTTTTAGTACTTGTTCTATTTGCATTGCGTAAGTAGGACATCCTCATATCCACCGGAACACACACACACATATATATATATATATATATATATACACACATATATACATAAATTCAATTTGGTGATCATTTTTTTGCATATATGGATGATCTTCATAAATGACCACACTCTTTATCAAAATTTACCCTACTATCTATGGTATCATTTTAGAAAATTTGAAATAGTATTTCAAAATATTGGATA
>chr1A:174688-175102
CTTCGTCGGTGCAGCCAACAGTCACTAAGTCCTTATTCAGCGAACGTAACCAATCATTGGCTACCATTGGGTCAGTTGTTCCAGAGAACCTTGCAGGGCTTAGTTTCAAAAACCTGGTGAGCATATCATGTGCTGGAGGATTATTCTGATTGTTGTTGTTGTTGTTCTGATTCTGATTTTGCAGGAACTGTGTCATCATCTGCATTAGTTGAGCTTGTCCTGTCAACATCTGCTGTATCTCGGGAGTGATGTTGTTGTTGTTGTTTTGGTTGGGGTTGGCACGACGGGGAGGCATATGTTGGTTTGATAGGGATAAGTAATGATAGATAGGCATAAGGAAGTTCTGACACACAACTCAAAACAACCAAGTTCAATCAAACAATCAATCAATCAACAAGCAACTATCCACTTAAC
>chr1A:182555-182866
```


Thiel T, Michalek W, Varshney RK, Graner A. Exploiting EST databases for the development and characterization of gene-derived SSR-markers in barley (Hordeum vulgare L.). Theor Appl Genet. 2003 Feb;106(3):411-22. doi: 10.1007/s00122-002-1031-0. Epub 2002 Sep 14. PMID: 12589540.
