
# Table of Contents
 <p><div class="lev1 toc-item"><a href="#Introduction" data-toc-modified-id="Introduction-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>Introduction</a></div><div class="lev2 toc-item"><a href="#FASTQ-files" data-toc-modified-id="FASTQ-files-11"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>FASTQ files</a></div><div class="lev2 toc-item"><a href="#Quality-control-of-the-raw-reads" data-toc-modified-id="Quality-control-of-the-raw-reads-12"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>Quality control of the raw reads</a></div><div class="lev2 toc-item"><a href="#Making-sense-of-reads" data-toc-modified-id="Making-sense-of-reads-13"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>Making sense of reads</a></div><div class="lev1 toc-item"><a href="#File-formats-explained" data-toc-modified-id="File-formats-explained-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>File formats explained</a></div><div class="lev2 toc-item"><a href="#BED-format" data-toc-modified-id="BED-format-21"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>BED format</a></div><div class="lev2 toc-item"><a href="#narrowPeak-format" data-toc-modified-id="narrowPeak-format-22"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>narrowPeak format</a></div><div class="lev2 toc-item"><a href="#bigWig" data-toc-modified-id="bigWig-23"><span class="toc-item-num">2.3&nbsp;&nbsp;</span>bigWig</a></div><div class="lev1 toc-item"><a href="#Entering-into-the-dark-side" data-toc-modified-id="Entering-into-the-dark-side-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Entering into the dark side</a></div><div class="lev2 toc-item"><a href="#Access-to-the-CRG-cluster" data-toc-modified-id="Access-to-the-CRG-cluster-31"><span class="toc-item-num">3.1&nbsp;&nbsp;</span>Access to the CRG cluster</a></div><div class="lev2 toc-item"><a href="#The-Unix-environment" data-toc-modified-id="The-Unix-environment-32"><span class="toc-item-num">3.2&nbsp;&nbsp;</span>The Unix environment</a></div><div class="lev2 toc-item"><a href="#Creating-a-directory" data-toc-modified-id="Creating-a-directory-33"><span class="toc-item-num">3.3&nbsp;&nbsp;</span>Creating a directory</a></div><div class="lev2 toc-item"><a href="#Getting-the-data" data-toc-modified-id="Getting-the-data-34"><span class="toc-item-num">3.4&nbsp;&nbsp;</span>Getting the data</a></div><div class="lev2 toc-item"><a href="#Moving-into-it" data-toc-modified-id="Moving-into-it-35"><span class="toc-item-num">3.5&nbsp;&nbsp;</span>Moving into it</a></div><div class="lev2 toc-item"><a href="#Listing-files-and-directories" data-toc-modified-id="Listing-files-and-directories-36"><span class="toc-item-num">3.6&nbsp;&nbsp;</span>Listing files and directories</a></div><div class="lev2 toc-item"><a href="#README" data-toc-modified-id="README-37"><span class="toc-item-num">3.7&nbsp;&nbsp;</span>README</a></div><div class="lev2 toc-item"><a href="#Viewing,-editing-and-executing-code" data-toc-modified-id="Viewing,-editing-and-executing-code-38"><span class="toc-item-num">3.8&nbsp;&nbsp;</span>Viewing, editing and executing code</a></div><div class="lev2 toc-item"><a href="#Relative-vs-absolute-paths" data-toc-modified-id="Relative-vs-absolute-paths-39"><span class="toc-item-num">3.9&nbsp;&nbsp;</span>Relative vs absolute paths</a></div><div class="lev1 toc-item"><a href="#Genome-arithmetic-I:-peaks-overlap" data-toc-modified-id="Genome-arithmetic-I:-peaks-overlap-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Genome arithmetic I: peaks overlap</a></div><div class="lev2 toc-item"><a href="#Objective" data-toc-modified-id="Objective-41"><span class="toc-item-num">4.1&nbsp;&nbsp;</span>Objective</a></div><div class="lev2 toc-item"><a href="#Data" data-toc-modified-id="Data-42"><span class="toc-item-num">4.2&nbsp;&nbsp;</span>Data</a></div><div class="lev2 toc-item"><a href="#Code" data-toc-modified-id="Code-43"><span class="toc-item-num">4.3&nbsp;&nbsp;</span>Code</a></div><div class="lev3 toc-item"><a href="#Difficult-bits-of-code-(optional)" data-toc-modified-id="Difficult-bits-of-code-(optional)-431"><span class="toc-item-num">4.3.1&nbsp;&nbsp;</span>Difficult bits of code (optional)</a></div><div class="lev2 toc-item"><a href="#Results" data-toc-modified-id="Results-44"><span class="toc-item-num">4.4&nbsp;&nbsp;</span>Results</a></div><div class="lev3 toc-item"><a href="#Inspecting-the-output" data-toc-modified-id="Inspecting-the-output-441"><span class="toc-item-num">4.4.1&nbsp;&nbsp;</span>Inspecting the output</a></div><div class="lev3 toc-item"><a href="#Counting-peaks" data-toc-modified-id="Counting-peaks-442"><span class="toc-item-num">4.4.2&nbsp;&nbsp;</span>Counting peaks</a></div><div class="lev3 toc-item"><a href="#Yes,-the-numbers-are-right..." data-toc-modified-id="Yes,-the-numbers-are-right...-443"><span class="toc-item-num">4.4.3&nbsp;&nbsp;</span>Yes, the numbers are right...</a></div><div class="lev1 toc-item"><a href="#Genome-arithmetic-II:-annotate-peaks" data-toc-modified-id="Genome-arithmetic-II:-annotate-peaks-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>Genome arithmetic II: annotate peaks</a></div><div class="lev2 toc-item"><a href="#Objective" data-toc-modified-id="Objective-51"><span class="toc-item-num">5.1&nbsp;&nbsp;</span>Objective</a></div><div class="lev2 toc-item"><a href="#Data" data-toc-modified-id="Data-52"><span class="toc-item-num">5.2&nbsp;&nbsp;</span>Data</a></div><div class="lev2 toc-item"><a href="#Code" data-toc-modified-id="Code-53"><span class="toc-item-num">5.3&nbsp;&nbsp;</span>Code</a></div><div class="lev2 toc-item"><a href="#Results" data-toc-modified-id="Results-54"><span class="toc-item-num">5.4&nbsp;&nbsp;</span>Results</a></div><div class="lev1 toc-item"><a href="#Acknowledgements" data-toc-modified-id="Acknowledgements-6"><span class="toc-item-num">6&nbsp;&nbsp;</span>Acknowledgements</a></div>

# Introduction

## FASTQ files

Simplistically speaking the output of high-throughput sequencing (HTS) machines are reads. These are stored in FASTQ files, which get this name from the [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format) in which reads are written. Essentially, a FASTQ file has millions of the following blocks:

**TODO: USE A FASTQ OF OURS**
```
@D00733:172:CA8BAANXX:8:2211:1214:2222 1:N:0:AACCAG
ATTGCCTAGGTTTTCTTCTAGGGTTTTTATGGTTTTAGGTCTAAAACACC
+
AB?ABFGFG>1?1FFGGFDGG1=1>/0C:111<EGGGG>1EFG@FF1FGF
```
Lines:
1. unique read ID plus sequencing information (Illumina machine, index...)
2. read sequence
3. irrelevant waste of space...
4. per-position quality


## Quality control of the raw reads

FASTQ files as obtained from the sequencer, typically referred to as "raw reads" are passed to programs which assess the quality of the reads in many ways. One of the most commonly used of such programs and the one we use is [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), which generates web-like reports (exampls of [good](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) and [bad](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html) data).


## Making sense of reads

Reads themselves do not tell much so there are some steps that are commonly performed in order to make sense of ChIP-seq data (many of them are common to many \*seq applications although each of them will require specific variations):
- trimming of technical sequences and low-quality ends
- alignment to the genome reference sequence
- generation of read counts profiles
- identification of regions of signal enrichment (aka peaks or binding sites)

These steps are illustrated [here](https://drive.google.com/open?id=0B-MXr-KyKmm6a1JOQnM0bUFSNGM).

# File formats explained

There are many HTS-related file formats out there (too many indeed!) so here we will only cover those used in this session.

## BED format

A tab-delimited no-header format to define the genomic coordinates and (optionally) some attributes of a set of genomic features. In its simplest form it has 3 mandatory fields:
```
chr22	1000	5000
chr22	2000	6000
```
Fields:
1. chromosome
2. feature start
2. feature end

A [BED file](https://genome.ucsc.edu/FAQ/FAQformat#format1) can have 9 additional fields/columns, which allow embedding additional attributes as we will see next.

## narrowPeak format

A format generated by the [ENCODE](https://www.encodeproject.org/) project to provide ChIP-seq peaks. It derives from the BED format and has the 3 first fields defining the genomic coordinates plus 7 more fields with information on the peaks (see [here](https://genome.ucsc.edu/FAQ/FAQformat#format12)). As we will see later, besides the genomic coordinates fields, the fields we pay more attention to are the 'signalValue' (signal enrichment) and the 'qValue' (significance).

## bigWig

More generally, *"the [bigWig format](https://genome.ucsc.edu/FAQ/FAQformat#format6.1) is for display of dense, continuous data"*. Here specifically is used to represent the number of read counts, preferably normalised, in genomic bins (i.e. ChIP-seq signal). Unlike the two previous formats, which can opened as a text file (yet this is not recommended as they can be very big!), bigWigs are binary files and cannot be opened.

# Entering into the dark side

## Access to the CRG cluster

Most of the analyses we will perform in this session can be run locally from your computer. However, there are several reasons to perform them in the CRG cluster:
- programs we will use (and other you may want to use in the future) are already installed
- more powerful
- lab data are stored there
- data are backed up

Open a Terminal (Applications > Utilities > Terminal): you are now in the command line! To connect to the CRG cluster by typing:
```
ssh ...
# password =
```

## The Unix environment

The Unix environment works very similar to the Windows/Mac/Ubuntu file systems:
- there are files and directories
- directories can have other files/directories within
- one can move down/up through directories
- files/directories can be created, copied, moved and deleted

We will use now Unix built-in functions/programs to both illustrate these concepts and prepare the data for this session. By the way, right now you are now in what it is called your *home* directory in the CRG cluster.


## Creating a directory

The built-in `mkdir` function will allow us to create a directory for this course:
```
mkdir bioinfo_course
```
Tips:
- avoid blank spaces; alternatives:
    - bioinfo_course
    - bioinfo-couser
    - bioinfo.course
    - bioinfoCourse
- avoid characters others than letters and numbers
- avoid file names starting with a number
- (personally) all lower-case


## Getting the data

The data for this session has been previously compiled in the `session01` directory, which is stored somewhere else in the CRG cluster. We will make a copy of the `session01` directory into the `bioinfo_course` you just created in your *home*.
```
cp -r /users/project/4DGenome/bioinfo_course/session01 bioinfo_course
```
- `cp` is the built-in function to copy files/directories
- `-r` is an option that modulates the behaviour of `cp` (options are explained in more detail later)
- `/users/project/4DGenome/bioinfo_course/session01`: path to directory that is being copied (it should be read as `session01` is within `bioinfo_course`, which is within `4DGenome`, which is within `project`, which is within `users`)
- `bioinfo_course`: directory in which `session01` will be copied


## Moving into it

We can move into the directory we just copied with `cd`:
```
cd bioinfo_course/session01
```


## Listing files and directories

Now can use the `ls` function to see the content of the current directory:
```
ls
```

The output of most Unix functions can be modulated by passing them different options. For instance, if we use:
```
ls -l -t -h
```
- `ls`: lists files and directories
- `l`: shows content as a list together with additional information (e.g. size, last time edited)
- `t`: sorts the date of the last edit
- `h`: shows size in human-readable units (i.e. KB, MB...)

The above can be compressed as:
```
ls -lth
```
You can learn about additional options with:
```
man ls
```
- man can be applied to any Unix function (we already saw it with `cp -v`)
- you can quit the manual page with by typing `q`

By default `ls` lists the content in the current directory, but it can also used to list the content in any directory:
```
ls -l data
```

And you can go as down as you want in the path, for instance:
```
ls -l analysis/overlap_peaks
```
This directory is still empty so nothing is shown.


## README

Inspecting the content within the `session01` directory reveals the `README.html` --the very same file you are reading now! We could have named this file with any name. However, including a [README](https://en.wikipedia.org/wiki/README) file is a convention computing people use to say "hey, here you will find information it is worth reading to start with". README files are normally just plain text (\*.txt) or [Markdown](https://guides.github.com/features/mastering-markdown/) (\*.md) files (the latter format enables some text formatting, which makes reading easier); here I have converted it to an HTML file for the sake of showing it online.


## Viewing, editing and executing code

So far we have executed code by either typing or pasting commands into the terminal. For convenience, code is typically packed into executable text files called scripts. The scripts used in this session are here:
```
ls -l scripts
```
The `*.sh` termination is a convention to indicate that they are unix shell executables files. We can see the code in the script with the `cat` function, which prints the content of a non-binary file:
```
cat scripts/hello_world.sh
```
- `#!/bin/bash`: this line is called the [shebang](https://en.wikipedia.org/wiki/Shebang_(Unix) and you will see it, or variations of it, at the beginning of in many scripts
- `echo "Hello, world"`: the `echo` function will print the quoted text that follows it ("Hello, world" in this case)

However, `cat` does not allow editing its content, something that can be done with a text editor. There are many text editors. I personally like [Sublime](https://www.sublimetext.com/) and I recommend you installing it for this session.

Finally, scripts can be executed just by typing its name:
```
hello_world.sh
```
You probably got the following error:
```
-bash: hello_world.sh: command not found
```
This is because `hello_world.sh` is not in your current directory but within `scripts`. We can execute it directly from the current directory, that is, no need to use `cd`, with:
```
scripts/hello_world.sh
```


## Relative vs absolute paths

Relative paths are paths **relative** to current directory, e.g.:
```
scripts/hello_world.sh
```
On the other hand, absolute paths include **all** the directories above (the `pwd` function shows the absolute path of the current directory); for instance:
```
/users/GR/mb/jquilez/bioinfo_course/session01/scripts/hello_world.sh
```


# Genome arithmetic I: peaks overlap


## Objective

The objective of this analysis is finding what peaks are shared by 2 ChIP-seq samples and what are sample-specific. 


## Data

The lists of peaks identified in each sample are in the narrowPeak format (see above) and can be found in:
```
ls -l data/*narrowPeak
```
We just introduced wilcards (`*`). In the command above we used `*` to list all files within `data` that ended with "narrowPeak", which, as you may have guessed, is the extension for the narrowPeak files.


## Code

The code needed to perform this analysis is in the following script:
```
ls -l scripts/overlap_peaks.sh
```
Let's go through the different steps performed by the script by looking at it. For the firs time I recommend you executing the code by copying small bits of code and pasting them onto the terminal. This is good to understand what each bit is doing. Alternatively, all the commands of a script are executed at a time by running the script:
```
scripts/overlap_peaks.sh
```

Regardless of how you execute the code, you should get that the number of filtered peaks is:
- sample 1: 191 peaks
- sample 2: 44,855 peaks


### Difficult bits of code (optional)

**Peaks filtering**

In the script above we used the following line of code to filter peaks based on their enrichment and significance:
```
awk -v min_qval=$min_qval -v min_enrichment=$min_enrichment '($7 > min_enrichment) && ($9 > min_qval)' $ifile1 | grep -v "chrM\|chrY\|chrUn\|luciferase" > $tfile1
```
- `awk` invokes the [awk programming language](https://en.wikipedia.org/wiki/AWK), which is very useful to deal with tabular files in Unix
- `-v min_qval=$min_qval`: `-v` is an awk option that allows defining awk-specific variables, which will be used later in the line:
    - the part that is to the left of the `=` is the awk-specific variable name
    - the part that is to the right of the `=` is the awk-specific variabla value, which can be a text string, a number or even the value of a Unix variable as we are using it here
    - variable name and value do not need to be the same...
- `-v min_enrichment=$min_enrichment`: another awk-specific variable is defined
- `'($7 > min_enrichment) && ($9 > min_qval)' $ifile1`: open `$ifile1` and print only rows that meet:
    - the value of column 7 > the `min_qval` awk-specific variable
    - the value of column 9 > the `min_enrichment` awk-specific variable
- `|`: this is called a pipe in Unix terms and it is a way to concatenate steps, in other words, here we are telling Unix "do not yet print the subset rows, we are gonna do something else"
- this *something else* is excluding rows containing some text
    - `grep` is another Unix built-in function to find rows that match a text/pattern
    - `-v` is a grep option that inverts the match, that is, find rows without a text/patern
    - `"chrM\|chrY\|chrUn\|luciferase"`: we can pass more than 1 text/patterns to match or avoid, in this case all the chromosomes we want to skip
- `>`: prints the rows left to the `$tfile1`

**Peaks overlap**

In the script above we used the following line of code to compare the genomic coordinates of the peaks in the 2 filtered peaks lists to find what peaks are shared and what peaks are specific of each list/sample:
```
bedtools multiinter -i $tfile1 $tfile2 -header -name sample1 sample2 > $ofile
```
- `bedtools` invokes the [bedtools](http://bedtools.readthedocs.io/en/latest/) toolset for genome arithmetic. bedtools not only is a very powerful tool in the analysis of genomic datasets but it is very well documented
- `multiinter`: bedtools comprises several tools and this is one of them, which finds whether regions in the files passed with `-i $tfile1 $tfile2` overlap or not
- `-header`: the output file will have a header (we will cover it later)
- `-name sample1 sample2`: arbitrary naming of some of the columns  of the header (we will cover it later)

## Results

### Inspecting the output

The output of the script is in the `analysis/overlap_peaks/overlap_peaks.bed` file. The `head` Unix built-in function allows inspecting the first 10 lines (by default, this can be changed with the `-n` option):
```
head analysis/overlap_peaks/overlap_peaks.bed
```
which prints on the screen:
```
chrom	start	end	num	list	sample1	sample2
chr1	818922	819164	1	sample2	0	1
chr1	877486	877706	1	sample2	0	1
chr1	906746	907037	1	sample2	0	1
chr1	1073709	1074296	1	sample2	0	1
chr1	1132612	1132766	1	sample2	0	1
chr1	1241579	1241766	1	sample2	0	1
chr1	1310495	1310844	1	sample2	0	1
chr1	1364222	1364570	1	sample2	0	1
chr1	1375146	1375481	1	sample2	0	1
```
- chrom: chromosome
- start: peak start
- end: peak end  
- num: number of samples in which the peak in the chrom:start-end coordinates is found
- list: samples in which the peak in the chrom:start-end coordinates is found
- sample1: whether the peak in the chrom:start-end coordinates is absent (0) or present (1) in sample1
- sample2: whether the peak in the chrom:start-end coordinates is absent (0) or present (1) in sample2


### Counting peaks

We can count now the number of peaks that are either shared or sample-specific:
```
grep -v chrom analysis/overlap_peaks/overlap_peaks.bed |cut -f5 |sort |uniq -c
```
Before we execute this command, let's break down its parts:
- as we are going to count rows in the `analysis/overlap_peaks/overlap_peaks.bed` file, we want to first exlcude the header; a way to do so is filtering out with `grep -v` that row containing `chrom`
- `cut -f5`: select values in the 5th column
- `sort`: sort them alphabetically
- `uniq -c`: identify unique values in the selected column and count how many times each of them is seen

By executing the command we can see that:
- shared peaks = 189
- sample1-specific peaks = 6
- sample2-specific peaks = 42,379
Here you have the numbers for a Venn diagram!

### Yes, the numbers are right...

If you do the maths the numbers do not make sense at first. Above we found that:
- sample 1: 191 peaks
- sample 2: 44,855 peaks

And:
- 189 + 6 is not 191
- 189 + 42,379 is not 44,855

This happens is due to (i) the way peaks coordinates are provided in the narrowPeak file produced by the peak caller MACS2 and (ii) how bedtools deals with it.

MACS2, after identifying peaks and their coordinates, can go for a second round with the [--call-summits option](https://github.com/taoliu/MACS) to *"deconvolve subpeaks within each peak called from general procedure. It's highly recommended to detect adjacent binding events. While used, the output subpeaks of a big peak region will have the same peak boundaries, and different scores and peak summit positions."*]. On the other hand, 2 peaks with exactly the same genomic coordinates will be considered as a single region by bedtools.

# Genome arithmetic II: annotate peaks


## Objective

The objective of this analysis is annotating a list of ChIP-seq peaks by finding the closest gene to each peak. 


## Data

The list of peaks identified in one of the samples:
```
ls -l data/gv_009_02_01_chipseq_peaks.narrowPeak
```

A BED file with the list of genes generated by the [GENCODE project](http://www.gencodegenes.org/) (GENCODE release 24):
```
ls -l data/gencode.v24.annotation_genes.bed 
```

With the following command we can see that `gencode.v24.annotation_genes.bed ` is a BED format with the 3 first mandatory columns specifying the genomic coordinates of the genes plus additional columns with other information: 
```
less -l data/gencode.v24.annotation_genes.bed 
```
`less` is another Unix built-in function which is very useful because:
- allows scrolling down the file (unlike `head`)
- combined with the `-S` option, the lines printed on the screen are not *broken* if their length exceeds the screen width

Fields:
1. chromosome
2. gene start
3. gene end
- gene ID
- *empty*
- strand
- source
- feature type (all are genes in this case)
- *empty*
- additional information separated by semicolons, e.g. gene_id, gene_type, gene_name, gene_status

## Code

The code needed to perform this analysis is in the following script:
```
ls -l scripts/annotate_peaks.sh
```
Let's go through the different steps performed by the script by looking at it. For the firs time I recommend you executing the code by copying small bits of code and pasting them onto the terminal. This is good to understand what each bit is doing. Alternatively, all the commands of a script are executed at a time by running the script:
```
scripts/annotate_peaks.sh
```

## Results

The output file should list all the filtered peaks along with the closest gene:
```
less -S analysis/annotate_peaks/gv_009_02_01_chipseq_peaks_annotated.bed
```
You will see a combination of colums from the peaks narrowPeaks file (left) plus columns from the genes files (right), obviously only for those peaks and genes which overlap.

From the previous analysis we know there are 191 filtered peaks in the sample we used now, so we can see whether this number of agrees with the number of rows of the output file:
```
wc -l analysis/annotate_peaks/gv_009_02_01_chipseq_peaks_annotated.bed
```
It shows there are 203 rows for 191 peaks. Can you guess why?

The answer is that in a few cases, 2 or more genes must be equally close to a given peak. We can find which peaks are present more than once with:
```
cut -f1-4 analysis/annotate_peaks/gv_009_02_01_chipseq_peaks_annotated.bed |uniq -d
```
- the code before the pipe (`|`) slices the first 4 columns of the file (rows with exactly the same values for these 4 columns will correspond to the same peak)
- the code after the pipe (`|`) identifies indeed which rows are repeated

There are 11 peaks linked to more than 1 gene:
```
chr1	44721703	44721849	gv_009_02_01_chipseq_peak_29
chr1	78004743	78004918	gv_009_02_01_chipseq_peak_35
chr11	62841402	62841684	gv_009_02_01_chipseq_peak_188
chr11	94708774	94708950	gv_009_02_01_chipseq_peak_214
chr12	100230022	100230180	gv_009_02_01_chipseq_peak_269
chr19	38831251	38831412	gv_009_02_01_chipseq_peak_577
chr19	40897286	40897497	gv_009_02_01_chipseq_peak_581
chr2	143335041	143335226	gv_009_02_01_chipseq_peak_635
chr22	29517780	29517921	gv_009_02_01_chipseq_peak_719
chr22	37887750	37887898	gv_009_02_01_chipseq_peak_726
chr7	139340600	139340853	gv_009_02_01_chipseq_peak_1053
```

Now we can look closer to the genes linked to the gv_009_02_01_chipseq_peak_29 peak:
```
grep -w gv_009_02_01_chipseq_peak_29 analysis/annotate_peaks/gv_009_02_01_chipseq_peaks_annotated.bed |less -S
```
Therefore, the peak overlaps with:
- the protein-coding gene *C1orf228*
- snRNA gene *RNU5F-1*

# Acknowledgements 
- [Enrique Vidal](https://github.com/qenvio) (Bioinformatician)
- [Gabriel Gonzalez](https://www.linkedin.com/in/gabriel-gonzalez-707758a?authType=NAME_SEARCH&authToken=HTr2&locale=en_US&trk=tyah&trkInfo=clickedVertical%3Amynetwork%2CclickedEntityId%3A31720183%2CauthType%3ANAME_SEARCH%2Cidx%3A1-1-1%2CtarId%3A1484914622561%2Ctas%3Agabriel%20gonz) (Bioinformatics Unit)
- [Damjana Kastelic](https://www.linkedin.com/in/damjana-kastelic-0668518b?authType=NAME_SEARCH&authToken=6otE&locale=en_US&trk=tyah&trkInfo=clickedVertical%3Amynetwork%2CclickedEntityId%3A321596430%2CauthType%3ANAME_SEARCH%2Cidx%3A1-1-1%2CtarId%3A1484914662749%2Ctas%3Adam) (Training Officer)
