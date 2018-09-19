# README #
This code computes Qabs and Qext from the refractive index or the dialectric constant using Mie's theory.        

Basic usage is
```python
from utils import get_q, load_eps
import matplotlib.pyplot as plt

# load refractive index data from file
data_eps = load_eps("eps_Sil.dat")

# compute qabs, qsca, qback for a grain of size asize
asize = 1e-7  # cm
data_q = get_q(data_eps, asize)

# plot results, data_q is a dictionary of wavelenght-dependent numpy arrays.
# note that wavelength is in micron
plt.loglog(data_q["wlen"], data_q["qabs"], label="$Q_{abs}$ this code")
plt.loglog(data_q["wlen"], data_q["qsca"], label="$Q_{sca}$ this code")

# add some ornaments to plot
plt.xlabel("$\\lambda / \\mu$m")
plt.ylabel("$Q_{sca}$ or $Q_{abs}$")
plt.legend(loc="best")
plt.savefig("output.png")

```