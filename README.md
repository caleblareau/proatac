<p align="left">
  <br><br><br>
  <img src="media/logo.png" width="50%"/>
</p>

[![Build Status](https://travis-ci.org/buenrostrolab/proatac.svg?branch=master)](https://travis-ci.org/buenrostrolab/proatac)
[![PyPI version](https://badge.fury.io/py/proatac.svg)](https://badge.fury.io/py/proatac)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

This package is maintained by [Caleb Lareau](mailto:caleblareau@g.harvard.edu) in the
[Buenrostro Lab](https://buenrostrolab.com). Source code is made freely available here
and a packaged install version is provided through [PyPi](https://pypi.python.org/pypi/proatac/).
<br><br>

## About<a name="about"></a>
The **proatac** package implements our data processing and quality control pipeline bulk
[ATAC-Seq](http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2688.html),
[single cell ATAC-Seq](http://www.nature.com/nature/journal/v523/n7561/full/nature14590.html),
and droplet-based ATAC-Seq data. 

## Table of Contents<a name="toc"></a>
- [About](#about)
- [Table of Contents](#toc)
- [Workflow Overview](#werk)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [.yaml Configuration](#yaml)

## Workflow Overview<a name="werk"></a>
NEED IMAGE HERE

## Installation<a name="installation"></a>
There are a few [dependencies](#dependencies) needed to get **proatac** to run. All are 
very common bioinformatics tools / languages and should be readily available in
most systems. However, **note that the current implementation of proatac is not supported
on Windows platforms**. 

Depending on your python environment, we generally recommend 

```
python3 -m venv venv3
source venv3/bin/active
pip3 install proatac
```

## Dependencies<a name="dependencies"></a>



- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [java](https://www3.ntu.edu.sg/home/ehchua/programming/howto/JDK_Howto.html)
- [macs2](https://github.com/taoliu/MACS)

We note that macs2 though also a PyPi package is only compatible with Python 2.7
whereas **proatac** is a Python 3 package. There's a good chance that macs2
is already living in your environment if you are reading this help page, which can
be tested using the following--

```
which macs2
```

and hopefully seeing a valid path. If not, one solution for macs2 install is to create
a separate python2 virtual environment using the following commands -- 

```
python2 -m venv venv
source venv/bin/active

pip install numpy
pip install wheel
pip install macs2
```

- [samtools](http://www.htslib.org/download/)

## .yaml configuration<a name="yaml"></a>



### Questions/comments/feedback
are always welcomed. Email [Caleb](mailto:caleblareau@g.harvard.edu) anytime! 
The easiest way for us to have correspondence (if appropriate/interesting
for the public) is through raising a [new issue](https://github.com/buenrostrolab/proatac/issues/new)
on the Github source. 


**proatac** logo made freely with [logomakr](https://logomakr.com/).
<br><br><br>