# README #
This code computes Qabs and Qext from the refractive index or the dialectric constant using Mie's theory.        

### Getting started ###
Basic usage is (see `test_01.py` file)
```python
from utils import Qabs_utils

# create object
q = Qabs_utils()

# load data from file with given format (here default format, but see next example)
q.load_eps("eps_Sil.dat")

# compute qabs and qsca for give grain size
asize = 1e-7  # cm
q.compute_q(asize)

# plot eps
q.plot(what=["real_eps", "im_eps"], fname="test_01_eps.png")

# plot refractive index
q.plot(what=["real_m", "im_m"], fname="test_01_m.png")

# plot qabs and qsca
q.plot(what=["qabs", "qsca"], fname="test_01_q.png")
```

### Different file formats ###
To load data from files with different formats, for example a file with wavelength, real and imaginary part of dielectric (see `test_03.py`):
```
from utils import Qabs_utils

# create object
q = Qabs_utils()

# load data from file with given format
q.load_eps("eps_CO.dat", labs=["wlen", "real_eps1", "im_eps"])
```
The file could be both space- or tab-separated.    

labs are     
`wlen`: wavelength in micron     
`real_eps1`: real part of dialectric - 1    
`real_eps`: real part of dialectric    
`im_eps`: imaginary part of dielectric     
`real_m1`: real part of refractive index - 1     
`real_m`: real part of refractive index     
`im_m`: imaginary part of refractive index

### Benchmark ###
The code has been benchmarked against the Original Astronomical Silicate (Draine & Lee 1984; Laor & Draine 1993) [link](https://www.astro.princeton.edu/~draine/dust/dust.diel.html).     
See `test_02.py`   
```
from utils import Qabs_utils

# create object
q = Qabs_utils()

# do benchmark
q.benchmark()
```


