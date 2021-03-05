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

base_dir = '/home/abzoghbi/data/mcg5-23-16/nustar_re_analysis'
```

    The autoreload extension is already loaded. To reload it, use:
      %reload_ext autoreload



```python
os.chdir(base_dir)
nu_data_dir = 'data/nustar'
os.system('mkdir -p %s'%nu_data_dir)
nustar_obsids = np.array(['60001046002', '60001046004', '60001046006', '60001046008'])
print('There are %d observations'%len(nustar_obsids))
print(', '.join(nustar_obsids))
```

    There are 4 observations
    60001046002, 60001046004, 60001046006, 60001046008



```python
os.chdir('%s/%s'%(base_dir, nu_data_dir))
az.data_tools.process_nustar_obsids(nustar_obsids)
```

    starting 60001046002 ...
    starting 60001046004 ...
    starting 60001046006 ...
    starting 60001046008 ...


## Spectral Extraction
-  Region size: 150’' for the source, and multiple regions for the background


```python
os.chdir('%s/%s'%(base_dir, nu_data_dir))
az.data_tools.extract_nustar_spec(nustar_obsids)
```

    starting 60001046002 ...
    starting 60001046004 ...
    starting 60001046006 ...
    starting 60001046008 ...



```python
# Summary of the spectra
nu_spec_data = az.data_tools.spec_summary(nustar_obsids, '%s_p/spec/spec_%d_a.grp')
```

    num   | obsid        | mjd_s      | mjd_e      | rat        | exp       
        1 | 60001046002  |  56446.354 |  56450.839 |       3.83 |     160.48
        2 | 60001046004  |  57068.498 |  57073.439 |       2.68 |      210.9
        3 | 60001046006  |   57074.01 |  57076.394 |       2.99 |     98.472
        4 | 60001046008  |  57094.779 |  57099.787 |       3.13 |     220.84



```python
# move the spectra to one place
os.system('mkdir -p data/spec')
for o in nustar_obsids:
    os.system('cp data/nustar/%s_p/spec/spec_* data/spec'%o)
```

## nustar light curves
### 1b: `3 79`


```python
os.chdir('%s/%s'%(base_dir, nu_data_dir))
lcdir, ebins, dt = '1b', '3 79', 512
az.data_tools.extract_nustar_lc(nustar_obsids, lcdir, ebins, dt)
```

### 1b: `3 10`


```python
os.chdir('%s/%s'%(base_dir, nu_data_dir))
lcdir, ebins, dt = '1c', '3 10', 512
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
