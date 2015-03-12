## Post-GWAS Prioritization through Integrated Analysis of Functional Annotation

### Dependencies
- [Python3.4](https://www.python.org/)
- [Requests](http://docs.python-requests.org/en/latest/)
- [Progress-Bar](https://pypi.python.org/pypi/progressbar2)
- [scipy](http://www.scipy.org)
- [numpy](http://www.numpy.org/)

### Usage:

```
prioritize3.py [-h] [-o DESTINATION_PATH] [-b NBINS] [-t THRESHOLD] [-a ANNOTATION_PATH] GWAS_DATA_PATH
```

See `prioritize3.py -h` for more detail

### DATA Format (for both GWAS_DATA and ANNOTATION):
A text file with n lines, each line contains chromosome number, coordinate and the GWAS p-value, separated by one tab (i.e. `'\t'`)

