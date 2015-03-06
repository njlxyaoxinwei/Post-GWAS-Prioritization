## Post-GWAS Prioritization through Integrated Analysis of Functional Annotation

### Dependencies
- [Python](https://www.python.org/)
- [Requests](http://docs.python-requests.org/en/latest/)
- [Progress-Bar](https://pypi.python.org/pypi/progressbar2)

### Usage:

```
prioritize3.py [-h] [-o DESTINATION_PATH] [-b NBINS] [-t THRESHOLD] GWAS_DATA_PATH
```

See `python prioritize3.py -h` for more detail

### GWAS_DATA Format:
A text file with n lines, each line contains chromosome number, coordinate and the GWAS p-value, separated by one tab (i.e. `'\t'`)

