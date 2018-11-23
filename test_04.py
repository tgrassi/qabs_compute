from utils import QabsManager
from optical import Optical
import numpy as np

# create object
q = QabsManager()


# load amorphous carbon, compute opacity. and plot
core_carbon = q.load_material("data/eps_carb_P93.dat", ["wlen", "real_m", "im_m"])
core_carbon.dust.rho_bulk = 2e0
core_carbon.compute_kappa()

# load and plot core
core_silicate = q.load_material("data/eps_Sil_Oss92.dat", ["wlen", "real_m", "im_m"])

# load and plot mantle
mantle_thin = q.load_material("data/eps_H93.dat", ["wlen", "real_m", "im_m"], units="1/cm")
vacuum = q.vacuum_as(mantle_thin)
# mantle_thin.extrapolate(7.99e2)
mantle_thin.add_impurity([core_carbon, vacuum], [0.0, 0.2]) #.11

mantle_thick = q.load_material("data/eps_H93.dat", ["wlen", "real_m", "im_m"], units="1/cm")
# mantle_thick.extrapolate(7.99e2)
mantle_thick.add_impurity([core_carbon, vacuum], [0.0, 0.2]) #.013

# mantle = q.load_material("data/eps_H93.dat", ["wlen", "real_m", "im_m"], units="1/cm")

#mantle.add_plot_ref_index("ref_index_mantle.png", ptype="plot", marker=".")

#mantle_thick.extrapolate(1e3)
#mantle.add_plot_ref_index("ref_index_mantle.png", ptype="plot", linestyle=":")


# compute core opacity and plot
q.clear_plots()

opt = Optical(None)

opt.load_kappa("data/kappa_oss94_v00.dat")
opt.add_plot_kappa("kappa.png", postfix=" (ref, ratio 0.0)", linestyle=":")

opt.load_kappa("data/kappa_oss94_v05.dat")
opt.add_plot_kappa("kappa.png", postfix=" (ref, ratio 0.5)", linestyle=":")

opt.load_kappa("data/kappa_oss94_v45.dat")
opt.add_plot_kappa("kappa.png", postfix=" (ref, ratio 4.5)", linestyle=":")

core_silicate.compute_kappa()


frac = np.array((1., 0.475862068965517))
frac /= sum(frac)

merged = q.merge_kappa([core_silicate, core_carbon], frac)
merged.add_plot_kappa("kappa.png", postfix=" (merged, ratio 0.0)",
                      linestyle="-")


# compute kappa for different shell thicknesses
for ii, vratio in enumerate([0.5, 4.5]):
    this_mantle = [mantle_thin, mantle_thick][ii]

    # make composite materials (Si+ice and aC+ice)
    composite_silicate = q.make_optical([core_silicate, this_mantle])
    composite_carbon = q.make_optical([core_carbon, this_mantle])

    aratio = ((vratio + 1e0)**(1./3.) - 1e0)
    # aratio = (vratio + 1e0)**(1./3.)
    composite_silicate.dust.aratio = aratio
    composite_silicate.compute_kappa()

    composite_carbon.dust.aratio = aratio
    composite_carbon.dust.rho_bulk = 2e0
    composite_carbon.compute_kappa()

    merged = q.merge_kappa([composite_silicate, composite_carbon], frac)
    merged.add_plot_kappa("kappa.png", postfix=" (merged, ratio %.1f)" % vratio,
                          linestyle="-") #, xlim=(1e2, 1e3), ylim=(1e0, 3e2))
