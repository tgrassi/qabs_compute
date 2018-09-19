# README #
This code computes Qabs and Qext from the refractive index or the dialectric constant using Mie's theory.        

Basic usage is (see `test_01.py` file)
```python
from utils import get_q, load_eps
import matplotlib.pyplot as plt

# load refractive index data from file
data_eps = load_eps("eps_Sil.dat")

# compute qabs, qsca, qback for a grain of size asize.
# note, asize in cm
asize = 1e-7  # cm
data_q = get_q(data_eps, asize)

# plot results, data_q is a dictionary of wavelenght-dependent numpy arrays.
# note that wavelength is in micron
plt.loglog(data_q["wlen"], data_q["qabs"], label="$Q_{abs}$")
plt.loglog(data_q["wlen"], data_q["qsca"], label="$Q_{sca}$")

# add some ornaments to plot
plt.xlabel("$\\lambda / \\mu$m")
plt.ylabel("$Q_{sca}$ or $Q_{abs}$")
plt.legend(loc="best")
plt.show()
```

To load data from files with different formats, for example a file with wavelength, real and imaginary part of dielectric:
```
from utils import get_q, load_eps
import matplotlib.pyplot as plt

# load refractive index data from file
data_eps = load_eps("your_file.dat", labs=["wlen", "real_eps", "im_eps"])
```
labs are     
wlen: wavelength in micron     
real_eps1: real part of dialectric - 1    
real_eps: real part of dialectric    
im_eps: imaginary part of dielectric     
real_m1: real part of refractive index - 1     
real_m: real part of refractive index     
im_m: imaginary part of refractive index
