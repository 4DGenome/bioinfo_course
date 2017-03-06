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
cd bioinfo_course/session02
mkdir -p analysis/{annotation,closest,overlap}
```

* `mkdir` means **m**a**k**e **dir**ectory
* `cd` **c**hange **d**irectory.
* The extra argument `-p` allows the system to create the **p**arent directories and do not complain if they already exist.
* `analysis/{annotation,closest,overlap}` is a compact notation that allows appending different elements (the ones inside braces) to a common prefix (`analysis/` in this case)

## Get the files

The data for this session is already in the cluster file system. You can copy to your working directory

```bash
cp -r /users/GR/mb/evidal/courses/bioinfo_course/session02/* .
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

There are several peak files corresponding to the experiments (`gv*_chipseq_peaks.bed`), one annotation bed file with information about genes (`gencode.v24.annotation_genes.bed`), the output tab-separated file of the differential expression analysis (`deseq2_untreated_0_progesterone_360.tsv`) and one file with the chromosome sizes (`chromosome.sizes`, in bp).

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
bedtools intersect -u -a $proms -b data/gv_009_01_01_chipseq_peaks.bed
```

This command subsets from the `$proms` file only those features overlapping peaks  in `data/gv_009_01_01_chipseq_peaks.bed`. We further use a pipe to count the number of promoters with PRBS

```bash
bedtools intersect -u -a $proms -b data/gv_009_01_01_chipseq_peaks.bed | wc -l
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
	basename $i .bed
	bedtools intersect -u -a $proms -b $i | wc -l
done
```

* `basename $i` takes one path and strips out all directories to get only the naked file name. It can take a second optional argument (`.bed` in this case) as a suffix, to remove it out also

Instead of just counting the number of promoters, we can store the resulting files

```bash
for i in data/gv*
do
	outfile=analysis/overlap/$(basename $i .bed).tss_overlap.bed
	bedtools intersect -u -a $proms -b $i > $outfile
done
```

* `outfile=analysis/overlap/$(basename $i .bed).tss_overlap.bed` is striping down the file name of `$i` (including the `.bed` suffix) and using it to create a new file name located in the `analysis/overlap/` sub-directory, storing it in the variable `outfile`.

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
bedtools shuffle -i data/gv_092_01_01_chipseq_peaks.bed -g $sizes
```
* `-i`, as usual, identifies the input bed file
* `-g`, as before, identifies the file with chromosome sizes

We can use the output of this command as the input peak file to the `intersect` command and check how many promoters have a fake-PRBS

```bash
bedtools shuffle -i data/gv_092_01_01_chipseq_peaks.bed -g $sizes | \
		 bedtools intersect -u -a $proms -b - | \
		 wc -l
```
and compare it with what we had before

```bash
wc -l analysis/overlap/gv_092_01_01_chipseq_peaks.tss_overlap.bed
```
It would be nice if we can automatize this for all peak files. Let's do it with a for loop!
 
```bash
for i in data/gv*
do
	name=$(basename $i .bed)
	tssoverlap=analysis/overlap/$name.tss_overlap.bed

	nreal=$(cat $tssoverlap | wc -l)
	nrandom=$(bedtools shuffle -i $i -g $sizes | \
	    bedtools intersect -u -a $proms -b - | wc -l)
					
	echo $name $nreal $nrandom

done
```

Wow! Maybe this code block is too much. Line by line:

* `name=$(basename $i .bed)`: "get the file name (without the sub-directories and the `.bed` suffix) and store it in the `name` variable
* `tssoverlap=analysis/overlap/$name.tss_overlap.bed`: create a file name from the `name` defined, adding the new sub-directories and the corresponding suffix (so, retrieve the original overlap file) and store it in the `tssoverlap` variable
* `nreal=$(cat $tssoverlap | wc -l)`: print the contents of `$tssoverlap` file, count the number of lines it has and store it in the variable `nreal`
* Pay attention to the definition of the `nrandom` variable: It is the very same command used in the previous for loop, but this time, instead of printing the output, we are telling the machine to store it in the `nrandom` variable
* Spit the original file name (`$name`), the real number of promoters with PRBS (`$nreal`) and the number of promoters with random fake-PRBS (`$nrandom`).


We can redirect the output to a results file (save it!).

```bash
for i in data/gv*
do
	name=$(basename $i .bed)
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
for i in data/gv*bed
do
	name=$(basename $i .bed)
	bedtools closest -d -t first -a $tss -b $i > analysis/closest/$name.closest_tss.bed
done
```
* `bedtools closest` takes two feature files (`-a` and `-b`) and returns, for each row in file `-a` the closest feature in file `-b`
* `-d` modifies the output so an extra column is added with the distance (in bp) to the corresponding closest feature.
* `t first` also modifies the output; it tells the tool to, whenever there is a tie (two PRBS at the same distance), just take the dirst one

# A "pinch" of R

We can use R to cross data sets, make some fancy plots and try summarize the information in some useful way. In order to launch the program, just type `R` :)

We will see a welcome message and that the prompt has changed (the prompt is the symbol at the very left of the line where we can type commands).

Once inside on `R` first thing we need to do is load some dependencies

```R
library("tidyverse")
options(stringsAsFactors = F)
```
Unlike `bash`, in R all functions/commands need to encapsulate the arguments between braces `(` `)`.

* `library("tidyverse")` allow us to load a set of handy functions that will ease data manipulation and ploting. There is a lot of information about them in the [web](http://tidyverse.org/)
* `options(stringsAsFactors = F)` sets the way R interpret files with letters (characters)

## Import data

The very first thing we need to do is import the data we want to work with. We will start setting the file names. We could do this "manually" coping the path to the file(s) of interest, but lazy guys hate typing long file names / paths, so they invented a some shortcuts to avoid this.

```R
difexp_file <- list.files("data",
                     pattern = "deseq2",
                     full = TRUE)
```
* The operator `<-` "assigns" the one value to one object. In this case, we are defining the object `difexp_file` as the output of the expression on the right hand side.
* `difexp_file` is just an arbitrary name, We could have name it differently, but never using a number as the first character.
* `list.files` is a function (note the braces) that takes as mandatory argument a path (i.e a directory) for which we want to **list** the **files** in.
* `pattern = "deseq2"` defines an optional argument that helps `list.files` to identify the files. In this case, we want to locate the files containing the word `"deseq2"`
* `full = TRUE` forces the function to return the full path of the files found, not just the file names.

So, in summary, the above command tells R to list all files in `"data"` directory that matches the pattern and output them as full paths. With this we store the file location of the results of the RNAseq differential expression analysis.

We can do something similar to locate the file containing the TSS

```R
tss_file <- list.files("analysis/annotation",
                       pattern = "tss",
                       full = TRUE)
```
Please note that this time we are looking in a different sub-directory.

Once R already knows which files to import, we have to read them.

```R
difexp <- read.delim(difexp_file)

tss <- read.delim(tss_file, header = FALSE)
names(tss) <- c("chr", "start", "end", "id", "score", "strand")
```

Again, we are creating new objects (`difexp` and `tss`) using the assign operator `<-`.

* `read.delim` is a function that reads a tab delimited file.
* The argument `header = FALSE` makes R to read a file without header (column names).

In the third line of this code chunk we are setting the attribute `names` of the object `tss`. This is needed, as the original file was headless.

* `c` is a function (note the braces) that concatenates a series of elements into one vector. Here we can see that words in R either have a special meaning (names assigned to object, functions or arguments) or have to be quoted.

## Describe data sets

The objects created with the `read.delim` function are `data.frames`, which is the way R stores data sets. A data frame is a collection of variables (encoded as vectors, in columns) for a series of registers / individuals (in rows).

We can inspect a `data.frame` in several ways:

* Showing the first rows

```R
head(difexp)
```
* Summarizing its columns

```R
summary(difexp)
```
* Getting its structure

```R
str(difexp)

str(tss)
```
In any case, we can check if the imported objects match our expectations.

## Basic interaction with

The `tidyverse` library provide us with some *intuitive* verbs to interact with our data. For instance, we can filter rows of a data frame

```R
difexp <- filter(difexp, !(is.na(padj)))
```
* `filter` function takes one `data.frame` as first argument, and then one (or many) vectors of conditions to filter for

Let's try to understand the condition from the inside to the outside of the braces: 

* `padj` is the name of a variable of the data set; it stores the adjusted p-value of the comparison
* `is.na` function evaluates each element of a vector (the `padj` variable) and returns a `TRUE` value if the element is `NA` (**N**ot **a**vailable, i.e. missing)
* the `!` operator switches logical values; it makes `TRUE` `FALSE` and vice versa

In a nutshell: only keep (`filter`) rows that do not present a missing value (`NA`) in the `padj` column

We can subset the columns of interest:

```R
tss <- select(tss, chr, id)
```
* `select` is a function that takes a `data.frame` and as many column names as we want to keep (`chr` and `id` here) and outputs just the selected columns.

We can check the result of these commands

```R
str(difexp)
str(tss)
```

## Merge information

Another interesting task that can be done is to "cross" information of two data sets, given that they have something in common.

In our example, `difexp` contains the differential expression results at gene level and `tss` contains the genomic location of the genes.

```R
str(difexp)
str(tss)

difexp <- inner_join(difexp, tss)

str(difexp)
```
* `inner_join` is a function that takes two `data.frame`s as arguments and returns another `data.frame` with all the columns of both original data sets and only the rows that they have in common

In order to work, both data sets must have some columns in common (in this case `id` variable, encoding the ENSEMBL gene IDs).

Please take a moment to compare the change in the structure of the object after including the genomic position information.

## Transform variables

Another useful "verb" is `mutate`. It is used to transform and create new variables inside a `data.frame` (first argument).

```R
difexp <- mutate(difexp,
                 direction = ifelse(log2FoldChange > 0, "up", "down"),
                 stat_diff = ifelse(padj < .01, "signf", "no-sig"),
                 bio_diff = ifelse(abs(log2FoldChange) > 1, "relevant", "irrelevant"),
                 diff_group = paste(stat_diff, bio_diff),
                 change = ifelse((stat_diff == "signf") & (bio_diff == "relevant"),
                                 "change", "no change"))

str(difexp)
```
In the above code we are creating four new variables for the data set `difexp`.

* `ifelse` is a function that evaluates its first argument (that should be a `TRUE` / `FALSE`) and returns, for each element of the evaluated vector, the second argument if the statement is `TRUE` and the third one else (if is `FALSE`)
* `abs` is a function that takes a numeric vector and returns its absolute value
* `paste` function creates a vector of characters out of two (or more) vectors 

For instance, the new variable `direction` will have the value `"up"` if `log2FoldChange > 0` and `"down"` otherwise. The variable `stat_diff` is `"signif"` if the adjusted p-values is below 0.01 (reflecting a statistical significant difference between conditions). The variable `bio_diff` is `"relevant"` if the relative difference, in absolute value, is above 2 times.  The variable `diff_group` just groups the possible outcomes of the previous `*_diff` variables. Finally, the variable `change` is "change" if and only if the difference is both statistically significant and biologically relevant.


## Ploting

We can create a volcano plot using `ggplot` functions

```R
ggplot(difexp,
       aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(alpha = .3)
```
* `ggplot` is a function than initializes a plot, usually taking two arguments: a `data.frame` and a aesthetics call
* `aes` is a function that allow us to tell to the `ggplot` function which variables we want to plot; in this case, we want to depict the log2 fold change on the X axis and minus the log10 of the adjusted p-value on the Y axis
* `geom_point` is one of the functions associated with `ggplot` that defines a layer to be represented in the plot; in this case, we want to depict some points
* `alpha` is a parameter that controls the transparency (0 fully transparent, 1 fully opaque)

`ggplot` is quite a powerful tool to make figures. For instance, we can color the plots according to the groups defined earlier just adding one extra argument to the `aes` call:

```R
ggplot(difexp,
       aes(x = log2FoldChange, y = -log10(padj),
           col = diff_group)) +
    geom_point(alpha = .3) +
    ylim(c(0, 50))
```
* `ylim` tells ggplot to restrict Y axis into the limits defined by the vector `c(0, 50)`

## Summarize

Having thousand of rows can be overwhelming. Why not try to summarize the information a little bit?

```R
summarize(difexp,
          n_up = sum(direction == "up"),
          n_tot = n(),
          p_up = n_up / n_tot,
          se_p_up = sqrt(p_up * (1 - p_up) / n()),
          ll = p_up - 2 * se_p_up,
          ul = p_up + 2 * se_p_up)
```

* `summarize` takes a data.frame and performs as many summaries of if as we can imagine.
* `n_up` is the sum of all rows where the column `direction` takes the value `"up"` (number of genes increasing their expression)

Can you guess what are the other columns?

We can combine several verbs in one expression to answer more complex questions. For instance, what is the proportion of genes increasing their expression value in chromosome chr13?

We could create a new object filtering the rows where `chr` takes the value `"chr13"` and then replicate the above summarizing command. Or we can take advantage of the pipe operator `%>%` that redirects the output of one function to the input of the next one. It is the `R` flavor of the UNIX `|` pipe.

```R
filter(difexp, chr == "chr13") %>%
    summarize(n_up = sum(direction == "up"),
              n_tot = n(),
              p_up = n_up / n_tot,
              se_p_up = sqrt(p_up * (1 - p_up) / n()),
              ll = p_up - 2 * se_p_up,
              ul = p_up + 2 * se_p_up)
```

## Grouping

Sometimes we want the results stratified. If we would like to know the proportion of genes going up for each chromosome, we could repeat last code batch changing the value of the chromosome, but it is more efficient (and more lazy) to do it all together.

```R 
group_by(difexp, chr) %>%
    summarize(n_up = sum(direction == "up"),
              n_tot = n(),
              p_up = n_up / n_tot,
              se_p_up = sqrt(p_up * (1 - p_up) / n()),
              ll = p_up - 2 * se_p_up,
              ul = p_up + 2 * se_p_up)
```

* `group_by` takes one `data.frame` and as many variables in it we want to stratify by
* `summarize` acting on a "grouped" `data.frame` returns the summaries per group

Maybe having the proportion of genes going up per chromosome is not very interesting. What about separate them whether the change is significant and relevant or not?

```R
group_by(difexp, change) %>%
    summarize(n_up = sum(direction == "up"),
              n_tot = n(),
              p_up = n_up / n_tot,
              se_p_up = sqrt(p_up * (1 - p_up) / n()),
              ll = p_up - 2 * se_p_up,
              ul = p_up + 2 * se_p_up)
```

So it seems the proportion of "real" changes in expression is smaller for up-regulated genes.

## Integrating more data

It could be interesting to add some information about PRBS.

We can select one of the TSS files annotated with the closest PRBS obtained before and import it into out session

```R
closest_files <- list.files("analysis/closest",
                        pattern = ".closest_tss.bed",
                        full = TRUE)
closest_files
closest_file <- closest_files[6]

```

This time there is more than one file matching the pattern. We can select an element out of a vector using the extractor operator (`[`) and the index of the element we want to extract.

* The `[6]` after the `closest_files` object indicates we just want the 6th element of the vector 

Once the file has been selected, we can proceed as before:

```R
closest <- read.delim(closest_file, header = F)
names(closest) <- c("chr_tss", "start_tss", "end_tss", "id",
                    "score_tss", "strand_tss",
                    "chr_peak", "start_peak", "end_peak", "dis")
```

And then create a new data set with both the differential expression analysis results and the proximity of PRBS to the TSS.

```R
dat <- inner_join(difexp, closest)

str(dat)
```

Let's explore the distribution of distances

```R
ggplot(dat,
       aes(x = dis / 1e3)) +
    geom_histogram(bins = 100) +
    xlim(c(0, 500)) +
	xlab("TSS distance to closest PRBS / Kbp")
```

* `geom_histogram` is another `ggplot` related functions; it creates an histogram with the variable assigned to the x axis in the `aes` call, with as many bins as defined in the `bins` argument
* `xlab` allow us to change the X axis label

## Relevant question?


> **Are those genes having a PRBS near their TSS more prone to suffer expression changes upon hormone treatment?** 

First, lets mark those genes with a peak overlapping the promoter and then we summarize

```R
mutate(dat,
       overlap = ifelse(dis < 5e3, "peak", "no peak")) %>%
    group_by(direction, overlap) %>%
    summarize(n_ch = sum(change == "change"),
              n_tot = n(),
              p_ch = n_ch / n_tot,
              se_p_ch = sqrt(p_ch * (1 - p_ch) / n()),
              ll = p_ch - 2 * se_p_ch,
              ul = p_ch + 2 * se_p_ch) %>%
    na.omit
```

What do you think?

### Bonus track (at your own risk)

We can go one step further and use some statistics to check if the probability of a gene to be regulated upon hormone treatment is associated with the presence of PRBS in the promoter and if this association cannot be explained by chance.

We will need to define an auxiliary matrix to get the 95 % CI out of the model

```R
alpha <- .05
wald_mat <- matrix(c(1, 0, 1, qnorm(alpha / 2), 1, qnorm(1 - alpha / 2)), 2)
colnames(wald_mat) <- c("est", "ll", "ul")
```
And then fit a logistic regression to model the probability for a gene to experiment a change (both significant and relevant) in terms of the direction of the change, the presence of PRBS in the promoter and their interaction

```R
mutate(dat,
       overlap = ifelse(dis < 5e3, 1, 0),
       down = ifelse(direction == "down", 1, 0)) %>%
    glm(change == "change" ~ down * overlap, binomial, .) %>%
    summary %$%
    coefficients %>%
    (function(x) x[,1:2] %*% wald_mat) %>%
    exp
```

Conclusions:

1. Overall, a gene has 80 % more chances to be down-regulated than up-regulated
2. A gene has 60 % more chances to be regulated if it has PRBS in its promoter
3. There is a non-additive effect of the direction and the presence of PRBS; a gene having a PRBS in its promoter has a probability 30 % higher to be down-regulated than what is expected considering both effects independently
 
# All together now

In case we want to replicate this steps without going through this whole document, the main commands for the `bash` part are gathered at [this script](scripts/lines.sh) and the `R` ones [here](scripts/lines.r).
