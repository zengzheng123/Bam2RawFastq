# Bam2RawFastq

##### Extract original FASTQ sequence for the alignments in BAM file

### Install:
Run "make" in the program directory to compile

### Usage:
```
[ARGUMENTS]

--bam         <string>      Input bam file
--fastq_list  <string>      Comma seperated list of raw fastq gzipped files. For pair end read, the list need to be in correct order: R1_file1,R2_file1,R1_file2, R2_file2,...
--output      <string>      Output file prefix, for PE read, the output files will be <prefix>_R1.fastq.gz & <prefix>_R2.fastq.gz. For single end read, the output file will be <prefix>.fastq.gz
--single_end                The script assumes the the data is pair end read by default, use this option if the data is single end read
--help                      Print command line usage

```

### FAQ:

```






``` 

### This software uses the following library

bamtools https://github.com/pezmaster31/bamtools

gzstream http://www.cs.unc.edu/Research/compgeom/gzstream/
