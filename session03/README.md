# Overview

# Access to the cluster

In [session01](https://github.com/4DGenome/bioinfo_course/tree/master/session01) and [session02](https://github.com/4DGenome/bioinfo_course/tree/master/session02) we worked from the command line in the [CRG cluster](http://www.linux.crg.es/index.php/Main_Page).

As a reminder, the CRG cluster is accessed with:
```bash
ssh -Y USER@ant-login.linux.crg.es
```
Where your USER should be provided by the system administrators of the CRG cluster and has the first letter of your name plus your (first) lastname, all lower-case (e.g. `jquilez` for Javier Quilez).

The command above brings us to our individual home directory in the cluster.

Most of the present session will be executed in [R](https://www.r-project.org/), which, as seen in [session02](https://github.com/4DGenome/bioinfo_course/tree/master/session02), can be run in the command line by typing:
```bash
R
```

[RStudio](https://www.rstudio.com/) offers a more user-friendly interface to R, and the CRG cluster hosts a RStudio server:
```
http://rstudio.linux.crg.es/
```
You will be asked for your USER (see above) and your password. By accessing to the RStudio server we will be working in a node of the CRG cluster from our home directory. You can get your home directory with:
```R
getwd()
```

<br>

# Load data

# Exploratory analysis

## Boxplots

##  PCA

# Differential expression

## P60 vs T0

## E60 vs T0

# Volcano plot

# Venn diagrams

# Metaprofiles
