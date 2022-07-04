# pynext

Python port of iNEXT.ind function from iNEXT 

**Code**: https://github.com/JohnsonHsieh/iNEXT/blob/master/R/iNEXT.r

**Paper**: https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12613

Descriptions: where count is a list of ints == counts/sample (counts is typically equivilant to "reads")
and endpoint is an integer meaning total number of counts or reads to extrapolate to. 
Uses a Baysian approach to interpolate and extrapolate species richness by applying Chao1 matrix calculations
Output is a dict of named lists, each contains n == # species objects for:

['reads', 'method', 'species', 'species_lower_conf',
       'species_upper_conf', 'coverage', 'coverage_lower_conf',
       'coverage_upper_conf']

**Usage**:

```
import pynext
import pandas as pd #optional 

count = [1, 5, 10, 17... ] 
endpoint = 500000

out = pynext(count, endpoint)
df  = pd.DataFrame
```
