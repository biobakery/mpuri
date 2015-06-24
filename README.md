# MPURI: Microbiome Plotting Utilities in R (version I) #

## Installation
Download the source code, either by using the links on the [downloads page](https://bitbucket.org/biobakery/mpuri/downloads) (try [v0.1 zip](https://bitbucket.org/biobakery/mpuri/get/v0.1.zip)) and unzip the file, or by cloning it using
 
```
#!

 hg clone ssh://hg@bitbucket.org/biobakery/mpuri
```

In R, run the following command to install dependencies:
```
#!r
install.packages(c("ggplot2", "stringr", "reshape2", "plyr"))
```

Then, run the following to install mpuri: 

```
#!r

install.packages("/path/to/mpuri", repos = NULL, type = "source")
```
Here, `/path/to/mpuri` is the location of the (unzipped) source code you downloaded previously. 

## Documentation

For documentation, please see the built-in R help system for mpuri.