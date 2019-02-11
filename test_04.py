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


# load mantle from Hudgins+1993
mantle_thin = q.load_material("data/eps_H93.dat", labels=["wlen", "real_m", "im_m"], units="1/cm")

# extrapolate mantle to wmax
# Note: the extrapolation algorithm is tuned for the current problem,
# and might get wrong results when applied to different data
mantle_thin.extrapolate(wmax)

# add spherical impurities to mantle_thin with core_carbon optical properties, and volume ratio 0.11
mantle_thin.add_impurity([core_carbon], [.11])


# load mantle to create another optical object
mantle_thick = q.load_material("data/eps_H93.dat", labels=["wlen", "real_m", "im_m"], units="1/cm")

# extrapolate (see comments above)
mantle_thick.extrapolate(wmax)

# add spherical impurities with core_carbon optical properties, and volume ratio 0.0113
mantle_thick.add_impurity([core_carbon], [0.013])


# set constant optical properties value from 5e1 micron to wmax micron
ice_catania = q.constant_material(5e1, wmax, (1.27 + 1j * 4e-3))

# load optical constant from material
ice_lab = q.load_material("data/eps_CO.dat", labels=["wlen", "real_eps", "im_eps"])


# clear plots
q.clear_plots()

# create optical object to load opacities from file (Ossenkopf+Henning 1994)
opt = Optical(None)

# load and plot bare grain opacity
opt.load_kappa("data/kappa_oss94_v00.dat")
opt.add_plot_kappa("kappa.png", postfix="OH94 ref, $v=0$", linestyle=":")

# load and plot 0.5 volume ratio opacity
opt.load_kappa("data/kappa_oss94_v05.dat")
opt.add_plot_kappa("kappa.png", postfix="OH94 ref, $v=0.5$", linestyle=":")

# load and plot 4.5 volume ratio opacity
opt.load_kappa("data/kappa_oss94_v45.dat")
opt.add_plot_kappa("kappa.png", postfix="OH94 ref, $v=4.5$", linestyle=":")


# compute bare silicate opacity
core_silicate.compute_kappa()

# compute silicate and carbon mixture ratio
frac = np.array((1., 0.475862068965517))
frac /= sum(frac)

# merge opacities and plot computed
merged = q.merge_kappa([core_silicate, core_carbon], frac)
merged.add_plot_kappa("kappa.png", postfix="OH94 comp, $v=0$",
                      linestyle="-", color="tab:blue")


# compute kappa for different shell thicknesses and materials
for ii, vratio in enumerate([0.5, 4.5, 0.5, 4.5, 0.5, 4.5]):
    # list of mantles optical properties
    this_mantle = [mantle_thin, mantle_thick, ice_catania, ice_catania, ice_lab, ice_lab][ii]

    # line properties
    linestyle = ["-", "-", "--", "--", "-", "-"][ii]
    color = ["tab:orange", "tab:green"][ii % 2]
    lab = ["OH94 comp", "OH94 comp", "BP98", "BP98", "this work", "this work"][ii]
    # label width
    if lab == "this work":
        lw = 3
    else:
        lw = 1

    # make composite materials (Si+ice and aC+ice)
    composite_silicate = q.make_optical([core_silicate, this_mantle])
    composite_carbon = q.make_optical([core_carbon, this_mantle])

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

    # add opacity to plot
    merged.add_plot_kappa("kappa.png", postfix="%s, $v=%.1f$" % (lab, vratio), color=color,
                          linestyle=linestyle, linewidth=lw)
