## Post-GWAS Prioritization through Integrated Analysis of Functional Annotation

### Dependencies
- [python version >=3.4 or >=2.7,6](https://www.python.org/)
- [requests](http://docs.python-requests.org/en/latest/)
- [progressbar2](https://pypi.python.org/pypi/progressbar2)
- [scipy](http://www.scipy.org)
- [numpy](http://www.numpy.org/)

### Python2 vs. Python3
prioritize2.py is the python2 compatible version, and prioritize3.py is the python3 version. 

### Usage

```
prioritize.py [-h] [-o DESTINATION_PATH] [-b NBINS] [-t THRESHOLD] [-a ANNOTATION_PATH] GWAS_DATA_PATH
```

See `prioritize.py -h` for more detail

### DATA Format
The following format is for both GWAS_DATA and ANNOTATION:

A text file with n lines, each line contains chromosome number, coordinate and the GWAS p-value, separated by one tab (i.e. `'\t'`)

### Build
freeze.py is used for building executables. Please use [cx_Freeze](http://cx-freeze.sourceforge.net/) for the build:

1. Make sure all dependencies are installed for the preferred python version with which you wish to run prioritize.

2. Rename the correct python source code to prioritize.py

3. Run (where `python` points to the preferred version of python):
```
python freeze.py build
```

The executable will be named `prioritize` under the `build` directory.
