# Introduction

## Overview of this session

The top rated topics in the mini-poll we did were:

* Genomic arithmetic
* Batch mode
* Advanced R

The first one was introduced in the [last session](https://github.com/4DGenome/bioinfo_course/blob/master/session01/README.md) and it will be expanded here.

The second one refers to automatizing routines, in order to avoid repetitive tasks, not only because of the boredom associated but also to avoid *operative errors* (human mistakes).

The third one comprises using `R` to manipulate data and try to get some sense out of it. As the **basic** `R` topic didn't get too much attention in the poll, we assumed certain degree of familiarity with the language ... just kidding, we will try to do *advanced* things but starting from the scratch (if that is possible).

## Enter the dark side

First thing we need to enter the dark side is to open a terminal. `Applications > Utilities > Terminal` (or something similar) would do the trick.

Then we need to access the cluster, just like we did [a couple of weeks ago](https://github.com/4DGenome/bioinfo_course/blob/master/session01/README.md#access-to-the-crg-cluster), but this time we will make some plots, so we have to tell the machine we want graphical forwarding with the argument `-Y`

```bash
ssh -Y USER@ant-login.linux.crg.es
```

Please insert your username instead of `USER`. The terminal will ask for your password. Type it and hit enter.

Once in the cluster the first thing we notice is a message (in red) warning us not to run heavy computation in the login nodes. The best way to go is request an interactive node in the queue system, where we can be as rude with the machine as we want to. This can be done with the simple command `qlogin`.

```bash
qlogin
```
However, depending of the status of the queue, the system may not provide us with an interactive mode immediately. During this session we are not going to run heavy routines, so we can work directly in the login node, thus avoiding to wait for the interactive node. **BUT** if you plan to play hard, please request one interactive node (or submit a job) so you don't get a complain from IT guys. 

In any case, please create a new directory for this session and change directory to it.

```bash
mkdir -p bioinfo_course/session02/
cp bioinfo_course/session02
mkdir -p analysis/{annotation,closest,overlap}
```

* `mkdir` means **m**a**k**e **dir**ectory
* `cd` **c**hange **d**irectory.
* The extra argument `-p` allows the system to create the **p**arent directories and do not complain if they already exist.
* `analysis/{annotation,closest,overlap}` is a compact notation that allows appending different elements (the ones inside braces) to a common prefix (`analysis/` in this case)

## Get the files

The data for this session is already in the cluster file system. You can copy to your working directory

```bash
cp -r /users/project/4DGenome/bioinfo_course/session02/* .
```

* `cp` stands for **c**o**p**y
* the extra argument `-r` allows to get **r**ecursively everything nested in the directory
* the wildcard `*` matches anything inside the source directory (including files and sub-directories)
* the `.` means "here" the current working directory and sets the destination of the files copied

## Configuration

As Xavi mentioned last week during his labmeeting, the operative system of the cluster has been updated recently. To make sure everything goes smoothly, please load the configuration file hosted by him.

```bash
source /users/GR/mb/jquilez/.bashrc
```
This will do some setup to make sure everything works on the new system.


# Description

The situation is similar to one in the previous session: We want to investigate the relationship between PRBS and genes in out favorite cell line (T47D). But this time we have a bunch of peak files coming from different experiments and we prefer not to copy the same command several times. Furthermore, we also have the results of a differential expression analysis comparing some RNAseq samples at T0 and after 6 hours of PR treatment.


We can have a look at the data, listing the contents of the corresponding directory

```bash
ls -lhr data
```

* `ls` stands for **l**i**s**t
* `-l` indicates **l**ong output
* `-h` indicates **h**uman readable file size units
* `-r` indicates **r**everse order

There are several peak files corresponding to the experiments (`gv*_chipseq_peaks.narrowPeak`), one annotation bed file with information about genes (`gencode.v24.annotation_genes.bed`), the output tab-separated file of the differential expression analysis (`deseq2_untreated_0_progesterone_360.tsv`) and one file with the chromosome sizes (`chromosome.sizes`, in bp).

You can have a look at any of this files using the `head` command, that print the first 10 lines of the file. For example

```bash
head data/gencode.v24.annotation_genes.bed
```

yields the first lines of the wide annotation file, containing some information or each gene (one per line) like the genomic location.

# Prepare annotation

Please have a look again at the gene annotation file. The genomic location refers to the start and end of the gene. For some analysis, it is interesting to focus on the transcription start sites (TSS) or the promoter regions.

## Set file paths as variables

To make the code more readable, easier to type thus less human-prone, it is a good practice to set some variables pointing to the files we are going to use.

```bash
genes=data/gencode.v24.annotation_genes.bed
sizes=data/chromosome.sizes
```

We can also use this to "future" files:

```bash
tss=analysis/annotation/gencode.v24.annotation_genes.tss.bed
proms=analysis/annotation/gencode.v24.annotation_genes.promoters.bed
```

## Get TSS

Obtaining the TSS is as simple as walking through the annotation file and storing the start position if the gene lies on the positive strand or the end position if it lies on the negative. This can be done in several ways, but a simple and fast one is using [`awk`](https://en.wikipedia.org/wiki/AWK) programming language.

```bash
awk -v OFS="\t" \
	'{if($6 == "+"){tss = $2}else{tss = $3}print $1, tss, tss, $4, $5, $6}' \
	$genes | \
	bedtools sort -i - > \
	$tss
```

Let's break the code in small pieces to understand it better. Firs of all, let's pay attention to some special characters at the end of the line:

* `\` at the end of the line is a escape character that allows us to break one single command into multiple lines. Here we used it only for aesthetic reasons. In fact, the machine doesn't care about the length of the command, but for humans is often easier to read a long command if it is broken down in several lines.
* `|` is the pipe operator in UNIX. Technically, it redirects the standard output of a command to the standard input of the following one, which means it recycles the output of one command as the input of the next one. It could seem strange at the beginning, but it is widely used as it saves a lot of computing time (there is no disk writing/reading, everything stays in the RAM)
* `>` is another special UNIX operator that redirects the standard output to a file; it writes the output to a file.

Please remember that we can obtain the value of a variable with the `$` operator. So `$genes` stands for the path to the file containing the gene annotation, as defined above.

The the `awk` code itself:

* `-v OFS="\t"` allows us to define a variable that can be used *inside* the `awk` call. In this case, we are setting the internal variable `OFS` (**o**utput **f**ield **s**eparator) to be `"\t"`, so we are telling `awk` to separate the columns of the output by tabs

The core of the `awk` call is surrounded by `'{` `}'`. In detail, we are telling the machine that if column number 6 (`$6` in `awk` code), corresponding to the strand, is equal to `"+"` then it should store the column number 2 (`$2` in `awk` code), corresponding to the start, in the variable `tss` and if not (`else`) it should store the column number 3 (`$3`), corresponding to the end, in the variable `tss`. It makes sense, doesn't it? Finally, we tell the program to print the values of the first column, then two times the `tss` value and then the columns 4, 5 and 6. The reason we print `tss` twice is because we want to obtain a [bed-like](https://github.com/4DGenome/bioinfo_course/blob/master/session01/README.md#bed-format) file (which has both start and end).

After that we make use of the sorting tool of [bedtools](http://bedtools.readthedocs.io/en/latest/), that sorts the bed file:
* `-i` flags the **i**nput file
* `-` indicates that the input argument comes from the standard input (the output of the previous command, through the pipe).

And voilà, we have a new file with the TSS of each gene. We can inspect it using head

```bash
head $tss
```

and compare it with the original annotation file to check that everything run as expected.

```bash
head $genes | cut -f 1-6
```

* the command `cut` "**cut**s" the input at the desired columns
* the argument `-f 1-6` indicates that we want the columns (**f**ields) form 1st to 6th

## Get promoters

Once we have the TSS, we can define the promoters as the 5 Kbp flanking regions of the TSS. Happily, bedtools has a [tool](http://bedtools.readthedocs.io/en/latest/content/tools/slop.html) that does the job.

```bash
bedtools slop \
		 -i $tss \
		 -g $sizes \
		 -b 5000 > \
		 $proms
```

* `slop` is the tool that increases the size of each feature in a bed file
* `-i` sets the **i**nput file
* `$tss` points to the TSS file produced before
* `-g $sizes` sets the **g**enomic size of chromosomes file (which is important as we don't want to extend a feature below 0 or above the total chromosome length)
* `-b 5000` tells bedtools to expand the feature 5000 bp on each side

# Overlap

##  Identify promoters with PRBS

We want to select those promoters that contain a PRBS.

In last session we used `bedtools multiinter` to identify common intervals among multiple interval files. Here we will follow a simpler approach using [`bedtools intersect`](http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html). It takes to `bed` files as arguments (`-a` and `-b`) and returns their intersection. The output can be modulated by several arguments. We can use the flag `-u` to output only those features in `-a` that overlap any feature of `-b`.

```bash
bedtools intersect -u -a $proms -b data/gv_009_01_01_chipseq_peaks.narrowPeak
```

This command subsets from the `$proms` file only those features overlapping peaks  in `data/gv_009_01_01_chipseq_peaks.narrowPeak`. We further use a pipe to count the number of promoters with PRBS

```bash
bedtools intersect -u -a $proms -b data/gv_009_01_01_chipseq_peaks.narrowPeak | wc -l
```

* `wc` is a command that performs **w**ord, line and byte **c**ounts of a file
* `-l` makes it only output the number of **l**ines

So far so good, but recall we have a bunch of peak files (`ls -lh data/gv*`) and we are lazy enough to type just the necessary. By the way, laziness is a desirable quality of a good (bio)informatician; it makes you automatize, saving time and possible errors. We can build a `for` loop that walks through all PRBS files, and for each of them, performs the overlap with promoters and outputs the number of promoters (an even save them into a file).

Baby steps first: Lets do a for loop that just outputs the file names

```bash
for i in data/gv*
do
	echo $i
done
```

* The first line begins the definition of the `for` statement, creating the variable `i` that will sequentially take the values contained in `data/gv*`, that is, all files in the sub-directory `data/` matching the pattern `gv*` (starting with `gv`),
* `do` in the second line tells the machine that the definition of the `for` statement has ended and now it starts a code block with the instructions that we want to execute on each iteration
* `echo $i` is the only instruction given here, asking to print the content of the variable `$i` (the peak file name).

As an exercise, include a instruction inside the `for` code block to count the lines of each file.

Now we can build a more complicated `for` loop to get the number of promoters that overlap a PRBS for each of the peak files. We are just recycling the `bedtools intersect` command used before, replacing the `-b` file by the file pointed by `$i`

```bash
for i in data/gv*
do
	basename $i
	bedtools intersect -u -a $proms -b $i | wc -l
done
```

* `basename $i` takes one path and strips all directories to get only the naked file name

Instead of just counting the number of promoters, we can store the resulting files

```bash
for i in data/gv*
do
	outfile=analysis/overlap/$(basename $i).tss_overlap.bed
	bedtools intersect -u -a $proms -b $i > $outfile
done
```

* `outfile=analysis/overlap/$(basename $i).tss_overlap.bed` is striping down the file name of `$i` and using it to create a new file name located in the `analysis/overlap/` sub-directory, storing it in the variable `outfile`.

We can encounter again the number of promoters bound by PRBS with a different instruction

```bash
wc -l analysis/overlap/*.tss_overlap.bed
```
and compare it with the original number of promoters

```bash
wc -l $proms
```
or the original number of peaks

```bash
wc -l data/gv*
```

## Are promoters enriched in PRBS?

Having the absolute numbers is nice, but it would be more useful to put this numbers in context. For instance, are PRBS prone to locate in promoters? More than what it would be expected by chance?

One way to answer this question is to create a fake random set of PRBS with the same characteristics (number and size) and check how many promoters have one of this fake-PRBS. Fortunately, `bedtools` has a tool to perform random shuffling of a given set of features(bed file).

```bash
bedtools shuffle -i data/gv_092_01_01_chipseq_peaks.narrowPeak -g $sizes
```
* `-i`, as usual, identifies the input bed file
* `-g`, as before, identifies the file with chromosome sizes

We can use the output of this command as the input peak file to the `intersect` command and check how many promoters have a fake-PRBS

```bash
bedtools shuffle -i data/gv_092_01_01_chipseq_peaks.narrowPeak -g $sizes | \
		 bedtools intersect -u -a $proms -b - | \
		 wc -l
```
and compare it with what we had before

```bash
wc -l analysis/overlap/gv_092_01_01_chipseq_peaks.narrowPeak.tss_overlap.bed
```
It would be nice if we can automatize this for all peak files. Let's do it with a for loop!
 
```bash
for i in data/gv*
do
	name=$(basename $i)
	tssoverlap=analysis/overlap/$name.tss_overlap.bed

	nreal=$(cat $tssoverlap | wc -l)
	nrandom=$(bedtools shuffle -i $i -g $sizes | \
	    bedtools intersect -u -a $proms -b - | wc -l)
					
	echo $name $nreal $nrandom

done
```

Wow! Maybe this code block is too much. Line by line:

* `name=$(basename $i)`: "get the file name (without the sub-directories) and store it in the `name` variable
* `tssoverlap=analysis/overlap/$name.tss_overlap.bed`: create a file name from the `name` defined, adding the new sub-directories and the corresponding suffix (so, retrieve the original overlap file) and store it in the `tssoverlap` variable
* `nreal=$(cat $tssoverlap | wc -l)`: print the contents of `$tssoverlap` file, count the number of lines it has and store it in the variable `nreal`
* Pay attention to the definition of the `nrandom` variable: It is the very same command used in the previous for loop, but this time, instead of printing the output, we are telling the machine to store it in the `nrandom` variable
* Spit the original file name (`$name`), the real number of promoters with PRBS (`$nreal`) and the number of promoters with random fake-PRBS (`$nrandom`).


We can redirect the output to a results file (save it!).

```bash
for i in data/gv*
do
	name=$(basename $i)
	tssoverlap=analysis/overlap/$name.tss_overlap.bed

	nreal=$(cat $tssoverlap | wc -l)
	nrandom=$(bedtools shuffle -i $i -g $sizes | \
	    bedtools intersect -u -a $proms -b - | wc -l)
					
	echo $name $nreal $nrandom

done > analysis/overlap/promoters_with_PRBD_and_random.txt
```

Finally, we can use some of the recently acquired `awk` magic to compute the enrichment

```bash
awk '{print $0, $2 / $3}' \
	analysis/overlap/promoters_with_PRBD_and_random.txt
```
`$0` stands for all columns, so we are telling `awk` to output all columns plus an extra one computed as the real over random number of promoters with PRBS.

*(Please take into account that the shuffling would yield different at your computer)*

So, are PRBS enriched in promoters???

# Closest PRBS to TSS

Now that you master the for-loop-way, let's use it to assign the closest PRBS to each TSS

```bash
for i in data/gv*
do
	name=$(basename $i)
	bedtools closest -d -a $tss -b $i > analysis/closest/$name.closest_tss.bed
done
```
* `bedtools closests` takes two feature files (`-a` and `-b`) and returns, for each row in file `-a` the closest feature in file `-b`
* `-d` modifies the output so an extra column is added with the distance (in bp) to the corresponding closest feature.
