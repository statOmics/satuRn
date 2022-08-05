# 1.4.1

We report a bug in satuRn 1.4.0. (Bioconductor release 3.15). The bug
was inadvertently introduced in satuRn 1.3.1 (from the former
Bioconductor devel). Note that the bug was not thus present in any of
the older Bioconductor releases 3.13 and 3.14 (satuRn 1.0.x, 1.1.x and
1.2.x).

Bug details:

Imagine a gene with three isoforms and two cell types. The goal is to assess DTU
between cell types. All isoforms are expressed in all cells of cell type 1. 
However, none of the isoforms are expressed in any of the cells in cell type 2
(i.e., the gene is not expressed in cell type 2). 

satuRn computes the log-odds of picking a certain isoform from the pool of 
isoforms in each cell type, and then compares these log-odds estimates between 
the cell types. However, in this example, the log-odds of picking a certain 
isoform from the pool of isoforms in cell type 2 cannot be computed, as there
is no data. Hence, the DTU test statistic should be *NA*. However, due to 
erroneous handling of *NA* estimates, which was inadvertently introduced in 
satuRn 1.3.1. while aiming to resolve github issue 16, the log-odds in cell 
type 1 will be compared to zero. Hence, (erroneous) results can be obtained for 
this contrast, even when there are no data in cell type 2.

Note that in many cases such isoforms may not pass filtering and would not get
evaluated altogether. However, when analyzing sprase scRNA-Seq datasets with
a lenient filtering criterium, this problem will apply, and will result in
mistakes in the inference.

# satuRn 1.4.0

satuRn version for Bioconductor release 3.15

# satuRn 1.3.1

* Bug fix: allow for fit errors to be propagated as NA results (github issue 15 by @jgilis)
* Bug fix: handle experimental designs with empty factor levels correctly (github issue 16 by @XueyiDong)
* Bug fix: identify transcripts that are the only expressed transcript of a gene and set NA results (github issue 17 by @jgilis)
* Bug fix: handle extreme z-scores in testDTU with diagplot2 option
* Enhancement:  plotDTU now allows for sparseMatrix input

# satuRn 1.1.2

* Improved error handling

# satuRn 1.1.1

* Initial Bioconductor release of satuRn


# satuRn 0.99.0

NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.

SIGNIFICANT USER-VISIBLE CHANGES

* Your main changes to a function `foo()` or parameter `param`.

BUG FIXES

* Your bug fixes. See more details at
<http://bioconductor.org/developers/package-guidelines/#news>.
