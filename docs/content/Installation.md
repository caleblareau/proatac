# Install stable version through PyPi
There are a few [dependencies](#dependencies) needed to get **proatac** to run. All are 
very common bioinformatics tools / languages and should be readily available in
most systems. However, **note that the current implementation of proatac is not supported
on Windows platforms**. 

Depending on your python environment, we generally recommend using a virtual environment
to keep python dependencies tidy. An example of installing **proatac** inside a new
python virtual environment called `venv3`

```
python3 -m venv venv3
source venv3/bin/active
pip3 install proatac
```

While installing **proatac** is obviously a great first step, make sure that all of the 
[ependencies](#dependencies) are met
