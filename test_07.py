# Use to reproduce the results 2026 paper
import matplotlib.pyplot as plt
from qabsmanager import QabsManager
from optical import Optical
import numpy as np
import matplotlib as mpl

plt.rcParams['text.usetex'] = True
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

q = QabsManager()


# load amorphous carbon optical data
core_carbon = q.load_material("data/eps_carb_P93.dat", labels=["wlen", "real_m", "im_m"])

# set carbon dust bulk density to 2 g/cm3, default is 2.9 g/cm3 (silicon)
core_carbon.dust.rho_bulk = 2e0

# compute opacity for bare carbon
core_carbon.compute_kappa()

# store the maximum wavelength found for carbon
wmax = max(core_carbon.materials[0].data["wlen"])

# load silicon material from file
core_silicate = q.load_material("data/eps_Sil_Oss92.dat", labels=["wlen", "real_m", "im_m"])


# load optical constant from material
#ice_lab_CO = q.load_material("data/eps_CO.dat", labels=["wlen", "real_eps", "im_eps"])

#ice_lab_CO2 = q.load_material("data/eps_CO2.dat", labels=["wlen", "real_eps", "im_eps"])


#fields = ["wlen", "real_eps", "im_eps", "dummy", "dummy"] ["wlen", "real_eps", "im_eps", "dummy", "dummy"]
fields = ["wlen", "dummy", "dummy", "real_m", "alpha"]

# f	eps_re	eps_im	n	alpha
ice_lab_8k = q.load_material("data/exp_2026_8K.dat", labels=fields, units="thz")
ice_lab_120k = q.load_material("data/exp_2026_120K.dat", labels=fields, units="thz")
ice_lab_150k = q.load_material("data/exp_2026_150K.dat", labels=fields, units="thz")

# clear plots
q.clear_plots()

# create optical object to load opacities from file (Ossenkopf+Henning 1994)
opt = Optical(None)

# load and plot bare grain opacity
opt.load_kappa("data/kappa_oss94_v00.dat")
opt.add_plot_kappa("kappa_oss94_v00.pdf", postfix="${\\rm OH94}$, $V$=$0$", linestyle="--",
                   linewidth=1.5, color="gray")

# load and plot 0.5 volume ratio opacity
#opt.load_kappa("data/kappa_oss94_v05.dat")
#opt.add_plot_kappa("kappa.pdf", postfix="OH94, $V=0.5$", linestyle=":",
#                   linewidth=2, color="tab:blue")

# load and plot 4.5 volume ratio opacity
opt.load_kappa("data/kappa_oss94_v45.dat")
#label = "${\\rm OH94}$ (${\\small \\rm H_2O}$:${\\small \\rm CH_3}$:${\\small \\rm OH}$:${\\small \\rm CO}$:${\\small \\rm NH_3}$\\,100:10:1:1) $V$=$4.5$"
label = "${\\rm OH94}$, $V$=$4.5$"
opt.add_plot_kappa("kappa_oss94_v45.pdf", postfix=label, linestyle=":",
                   linewidth=1.5, color="gray")

# compute bare silicate opacity
core_silicate.compute_kappa()

# compute silicate and carbon mixture ratio
frac = np.array((1., 0.475862068965517))
frac /= sum(frac)

# compute kappa for different shell thicknesses and materials
for ii, vratio in enumerate([4.5, 4.5, 4.5]):
    # list of mantles optical properties
    ice_mantle = [ice_lab_8k, ice_lab_120k, ice_lab_150k][ii]

    # line properties
    linestyle = ["-", "-", "-"][ii]
    # get colormap viridis
    lab = ["This work, 8K", "This work, 120K", "This work, 150K"][ii]

    lab = "${\\rm " + lab.replace(" ", "\\ ") + "}$"

    # make composite materials (Si+ice and aC+ice)
    composite_silicate = q.make_optical([core_silicate, ice_mantle])
    composite_carbon = q.make_optical([core_carbon, ice_mantle])


    # define size ratio
    aratio = ((vratio + 1e0)**(1./3.) - 1e0)

    # set size ratio silicon/ice material
    composite_silicate.dust.aratio = aratio
    composite_silicate.compute_kappa()

    # set size ratio carbon/ice material
    composite_carbon.dust.aratio = aratio
    composite_carbon.dust.rho_bulk = 2e0  # g/cm3
    composite_carbon.compute_kappa()


    # merge carbon and silicate opacities
    merged = q.merge_kappa([composite_silicate, composite_carbon], frac)
    merged.save_kappa("kappa_%s_%.1f.dat" % (lab, vratio))

    # add opacity to plot
    merged.add_plot_kappa("kappa.pdf", postfix="%s, $V$=$%.1f$" % (lab, vratio), color=None,
                          linestyle=linestyle, linewidth=2, xlim=(1e2, 999.), ylim=(5e-1, 2e2))

from matplotlib.ticker import FuncFormatter, ScalarFormatter


clight = 2.99792458e10  # speed of light in cm/s


def myFormatter(x, pos):
    return "%.0f" % x

def myFormatter2(x, pos):
    if abs(x - 1) < 1e-5:
        x = 1
    return "%.1f" % x

mfmt = FuncFormatter(myFormatter)
mfmt2 = FuncFormatter(myFormatter2)

ax = plt.gca()
#ax.xaxis.set_major_formatter(ScalarFormatter())
ax.xaxis.set_major_formatter(mfmt)
ax.xaxis.set_minor_formatter(mfmt)

ax = plt.gca().secondary_xaxis('top', functions=(lambda x:clight / (x * 1e-4) / 1e12, lambda x:1e12 * (x * 1e-4) / clight))
ax.set_xlabel("$\\nu$, THz")
ax.xaxis.set_major_formatter(mfmt2)
ax.xaxis.set_minor_formatter(mfmt2)
plt.tight_layout()
plt.savefig("kappa.pdf")
plt.show()