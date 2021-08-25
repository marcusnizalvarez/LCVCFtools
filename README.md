# LCVCFtools
LCVCFtools is a simple C++ program designed for working with VCF v4.2 files generated from low-coverage whole genome sequencing. This tools is useful for data optimization before genotype imputation steps by removing samples or variantes with excessive missing data.

# Compiling
You can easily compile this program using LCVCFtools.pro file and qmake tool, just run the following command:
```
sudo apt install qt5-qmake libboost-all-dev libz-dev  
git clone https://github.com/marcusnizalvarez/LCVCFtools.git  
cd LCVCFtools/  
qmake && make
```
# Usage
## Input mode 
* --vcf <STRING>           Read from VCF file. Use - to read from stdin.
* --gzvcf <STRING>         Read from Gzip compressed VCF file. Use - to read from stdin.
 
## Filter parameters
* --minGQ <INT>            Minimum genotype quality in PhredScale. [Default=20]
* --minDP <INT>            Minimum depth. [Default=5]
* --MAF <FLOAT>            Minor allele frequency, based on allele depth (AD). [Default=0.1]
* --minGCR <FLOAT>         Minimum genotype call rate [Default=0].
* --minDPR <INT> <FLOAT>   Minimum DP rate. Can be defined multiple times.
* --minGQR <INT> <FLOAT>   Minimum GQ rate. Can be defined multiple times.
 
 ## Other arguments
* --remove <STRING>        Remove samples from a list file (One sample name per line).
* --sample-stats           Output sample statistics to 'stats.tsv'.
* --keep-multiallelic      Don't skip multiallelic variants.
* --ID                     Generate generic ID, useful for programs like Plink.
* --verbose                Verbose mode.
* --help                   Print this message.
  
## Usage example  
```
./LCVCFtools --gzvcf example.vcf.gz - --minGQ 20 --minDP 5 --minGCR 0.25 --minDPR 5 0.5 --MAF 0.1 --sample-stats | gzip -c > output.vcf.gz
```
# Credits
Author: Marcus Vinicius Niz Alvarez  
Email: marcus.alvarez@unesp.br  
São Paulo State University, UNESP - Biotechnology Institute and Bioscience Institute, Botucatu, 18618-689, Brazil.

# License
LCVCFtools is under GNU GPLv3.0 license and Boost Library is under Boost Software License v1.0.
