# Install stable version through PyPi
There are a few [dependencies](http://proatac.readthedocs.io/en/latest/content/Dependencies.html)
needed to get **proatac** to run. All are 
very common bioinformatics tools / languages and should be readily available in
most systems. However, **note that the current implementation of proatac is not supported
on Windows platforms**. 

Depending on your python environment, we generally recommend using a virtual environment
to keep python dependencies tidy. An example of installing **proatac** inside a new
python virtual environment called `venv3` using the following sequence of commands--

```
python3 -m venv venv3
source venv3/bin/activate
```

The installation can then be specified using the following:

```
pip3 install proatac
```

# Install via GitHub

Though **not recommended**, a bleeding-edge (development) version can be installed
directly from Git. Again using a virtual environment--

```
python3 -m venv venv3
source venv3/bin/activate
pip3 install git+ssh://git@github.com/buenrostrolab/search/tree/master/proatac
```

While installing **proatac** is obviously a great first step, make sure that all of the 
[dependencies](http://proatac.readthedocs.io/en/latest/content/Dependencies.html) are met. 
Check out the next page for more detail. 

# Installing proatac as an environment module.

A common use case of **proatac** will be processing ATAC-seq data in a high-performance computing cluster environment.
As each computing cluster is different, you're probably best off inquiring with you sysadmin how to install **proatac**
Here are a few general tips though (modified from the [MultiQC installation guide](http://multiqc.info/docs/#installing-as-an-environment-module)):

A typical installation procedure with an environment module Python install might look
like this: (Note that `$PYTHONPATH` must be defined before pip installation; this can be 
specified by creating a virtual environment and loading it as shown above)

```
# Create the proatac version (replace 0.3 with whatever the current version is)
VERSION=0.3
INST=/path/to/software/proatac/$VERSION
module load python/3.6.1
mkdir $INST
export PYTHONPATH=$INST/lib/python3.6/site-packages
pip3 install --install-option="--prefix=$INST" multiqc
```

Once **proatac** is installed, we need to create a new module file.
While these vary between systems, but here's an example that also imports the necessary 
dependencies:

```
#%Module1.0#####################################################################
##
## proatac
##

set components [ file split [ module-info name ] ]
set version [ lindex $components 1 ]
set modroot /path/to/software/proatac/$version

proc ModulesHelp { } {
    global version modroot
    puts stderr "\proatac - use proatac $version"
    puts stderr "\n\tVersion $version\n"
}
module-whatis   "Loads proatac environment."

# load required modules
module load python/3.6.1
module load samtools/1.5.0
module load R/3.4.0
module load java/1.6.0
module load macs2/2.1.1.20160309
module load bowtie2/2.3.1

# only one version at a time
conflict proatac

# Make the directories available
prepend-path    PYTHONPATH  $modroot/lib/python3.6/site-packages
```


