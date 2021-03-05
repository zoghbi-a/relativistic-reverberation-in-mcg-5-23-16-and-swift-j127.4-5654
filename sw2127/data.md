Using: 
- ``heasoft 6.27.2`` 
- ``xspec 12.11.0k`` 
- NuSTAR caldb: 20200912


```python
import numpy as np
import os
import aztools as az
import matplotlib.pyplot as plt

%load_ext autoreload
%autoreload 2

base_dir = '/home/abzoghbi/data/swift_j2127.4p5654/nustar_re_analysis'
```


```python
os.chdir(base_dir)
nu_data_dir = 'data/nustar'
os.system('mkdir -p %s'%nu_data_dir)
nustar_obsids = np.array(['60001110002', '60001110003', '60001110005', '60001110007', '60402008002', 
                          '60402008004', '60402008006', '60402008008', '60402008010'])
print('There are %d observations'%len(nustar_obsids))
print(', '.join(nustar_obsids))
```

    There are 9 observations
    60001110002, 60001110003, 60001110005, 60001110007, 60402008002, 60402008004, 60402008006, 60402008008, 60402008010



```python
os.chdir('%s/%s'%(base_dir, nu_data_dir))
az.data_tools.process_nustar_obsids(nustar_obsids)
```

    starting 60001110002 ...
    starting 60001110003 ...
    starting 60001110005 ...
    starting 60001110007 ...
    starting 60402008002 ...
    starting 60402008004 ...
    starting 60402008006 ...
    starting 60402008008 ...
    starting 60402008010 ...


## Spectral Extraction
-  Region size: 150’' for the source, and multiple regions for the background


```python
os.chdir('%s/%s'%(base_dir, nu_data_dir))
az.data_tools.extract_nustar_spec(nustar_obsids)
```

    starting 60001110002 ...
    starting 60001110003 ...
    starting 60001110005 ...
    starting 60001110007 ...
    starting 60402008002 ...
    starting 60402008004 ...
    starting 60402008006 ...
    starting 60402008008 ...
    starting 60402008010 ...



```python
# Summary of the spectra
nu_spec_data = az.data_tools.spec_summary(nustar_obsids, '%s_p/spec/spec_%d_a.grp')
```

    num   | obsid        | mjd_s      | mjd_e      | rat        | exp       
        1 | 60001110002  |  56235.737 |  56236.754 |      0.836 |     49.202
        2 | 60001110003  |  56236.755 |  56237.327 |       1.06 |     28.765
        3 | 60001110005  |  56237.758 |   56239.28 |       1.14 |     74.583
        4 | 60001110007  |  56239.711 |  56240.561 |       1.25 |      42.11
        5 | 60402008002  |  58304.691 |  58306.338 |       1.96 |     71.255
        6 | 60402008004  |  58315.297 |  58316.814 |       1.46 |     71.568
        7 | 60402008006  |  58329.864 |  58331.508 |       1.38 |      72.13
        8 | 60402008008  |  58375.703 |  58377.285 |      0.967 |     72.875
        9 | 60402008010  |  58482.343 |  58483.861 |       1.17 |     74.246


## nustar light curves
### 1b: `3 79`


```python
os.chdir('%s/%s'%(base_dir, nu_data_dir))
lcdir, ebins, dt = '1b', '3 79', 512
az.data_tools.extract_nustar_lc(nustar_obsids, lcdir, ebins, dt)
```

### 22l3: two sets of log space between 3-10; 10-80
`np.concatenate([np.logspace(np.log10(3), 1, 14), np.logspace(1, np.log10(79), 14)[1:]])`

The last four bins are result of merging every two neighbouring bins (e.g. 8 of them)


```python
os.chdir('%s/%s'%(base_dir, nu_data_dir))
lcdir, ebins, dt = '22l3', ('3 3.3 3.6 4 4.4 4.8 5.2 5.7 6.3 6.9 7.6 8.3 9.1 10 '
                           '11.7 13.8 16.2 19 22 31 42 58 79'), 512
az.data_tools.extract_nustar_lc(nustar_obsids, lcdir, ebins, dt)
```


```python

```
