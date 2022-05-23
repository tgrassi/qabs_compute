# Use to reproduce the results from Gavdush+2022

from qabsmanager import QabsManager
from optical import Optical
import numpy as np

# create manager object
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
ice_lab_CO = q.load_material("data/eps_CO.dat", labels=["wlen", "real_eps", "im_eps"])

ice_lab_CO2 = q.load_material("data/eps_CO2.dat", labels=["wlen", "real_eps", "im_eps"])

# clear plots
q.clear_plots()

# create optical object to load opacities from file (Ossenkopf+Henning 1994)
opt = Optical(None)

# load and plot bare grain opacity
opt.load_kappa("data/kappa_oss94_v00.dat")
opt.add_plot_kappa("kappa.pdf", postfix="OH94, $V=0$", linestyle=":",
                   linewidth=2, color="k")

# load and plot 0.5 volume ratio opacity
#opt.load_kappa("data/kappa_oss94_v05.dat")
#opt.add_plot_kappa("kappa.pdf", postfix="OH94, $V=0.5$", linestyle=":",
#                   linewidth=2, color="tab:blue")

# load and plot 4.5 volume ratio opacity
opt.load_kappa("data/kappa_oss94_v45.dat")
opt.add_plot_kappa("kappa.pdf", postfix="OH94, $V=4.5$", linestyle=":",
                   linewidth=2, color="tab:orange")

# compute bare silicate opacity
core_silicate.compute_kappa()

# compute silicate and carbon mixture ratio
frac = np.array((1., 0.475862068965517))
frac /= sum(frac)

# compute kappa for different shell thicknesses and materials
for ii, vratio in enumerate([4.5, 4.5]):
    # list of mantles optical properties
    ice_mantle = [ice_lab_CO, ice_lab_CO2][ii]


    # line properties
    linestyle = ["-", "--"][ii]
    color = ["tab:orange", "tab:orange"][ii % 2]
    lab = ["This work, CO", "This work, CO$_2$"][ii]

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
    merged.add_plot_kappa("kappa.pdf", postfix="%s, $V=%.1f$" % (lab, vratio), color=color,
                          linestyle=linestyle, linewidth=1.5, xlim=(3e1, 1.1e3), ylim=(5e-1, 2.5e3))
