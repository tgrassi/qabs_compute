# README #
This code computes Qabs, Qsca, and Qext from the refractive index or the dialectric constant using Mie's theory.     
It also compute opacity for coated materials.      

### Getting started ###
Basic usage is (see `test_01.py` file)
```python
from qabsmanager import QabsManager

# create object
q = QabsManager()


# load core data from file
core = q.load_material("data/eps_carb_P93.dat", labels=["wlen", "real_m", "im_m"])

# compute opacity for core material
core.compute_kappa()

# print some info to screen
q.report()

# save opacity to plot
core.plot_kappa("kappa_01.png")

```

### Different file formats ###
To load data from files with different formats, for example a file with wavelength, real and imaginary part of dielectric:
```
from qabsmanager import Qabs_utils

# create object
q = Qabs_utils()

# load data from file with given format
q.load_eps("eps_CO.dat", labels=["wlen", "real_eps1", "im_eps"])
```
The file can be both space- or tab-separated.    

labs are     
`wlen`: wavelength in micron     
`real_eps1`: real part of dialectric - 1    
`real_eps`: real part of dialectric    
`im_eps`: imaginary part of dielectric     
`real_m1`: real part of refractive index - 1     
`real_m`: real part of refractive index     
`im_m`: imaginary part of refractive index

### Composite material ###
It is possible to use coated materials by loading their optical properties (see `test03.py`).
```python
from qabsmanager import QabsManager

# create object
q = QabsManager()

q.clear_plots()

# load core data from file
core = q.load_material("data/eps_carb_P93.dat", labels=["wlen", "real_m", "im_m"])

# load ice data from file
mantle = q.load_material("data/eps_H93.dat", labels=["wlen", "real_m", "im_m"])

# create new material from core and mantle
composite = q.make_optical([core, mantle])

# set size ratio silicon/ice material
composite.dust.aratio = 0.15

# compute opacity of the composite material
composite.compute_kappa()

# plot opacity
composite.add_plot_kappa("kappa_03.png", postfix="coated", xlim=(5e1, 1e3), ylim=(1, 1e3))

```


### Benchmark ###
The code has been benchmarked against the Original Astronomical Silicate (Draine & Lee 1984; Laor & Draine 1993) [link](https://www.astro.princeton.edu/~draine/dust/dust.diel.html).     
See `test_02.py`   
```
from qabsmanager import Qabs_utils

# create object
q = Qabs_utils()

# do benchmark
q.benchmark()
```


