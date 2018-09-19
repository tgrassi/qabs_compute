from utils import get_q, load_qref, load_eps
import matplotlib.pyplot as plt

# load refractive index data from file
data_eps = load_eps("eps_Sil.dat")

# compute qabs, qsca, qback for a grain of size asize
asize = 1e-7  # cm
data_q = get_q(data_eps, asize)

# load qabs data for a 1e-3 micron grain from file for comparison
data_qref = load_qref("Sil_21_1e3.dat")

# plot data to compare qabs
plt.loglog(data_q["wlen"], data_q["qabs"], label="$Q_{abs}$ this code")
plt.loglog(data_qref["wlen"], data_qref["qabs"], ":", label="$Q_{abs}$ Draine")

# plot data to compare qsca
plt.loglog(data_q["wlen"], data_q["qsca"], label="$Q_{sca}$ this code")
plt.loglog(data_qref["wlen"], data_qref["qsca"], ":", label="$Q_{sca}$ Draine")

# add some ornaments to plot
plt.xlabel("$\\lambda / \\mu$m")
plt.ylabel("$Q_{sca}$ or $Q_{abs}$")
plt.legend(loc="best")
plt.savefig("test_01.png")


