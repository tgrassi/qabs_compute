from utils import QabsManager
from optical import Optical
import numpy as np

# create object
q = QabsManager()

# load and plot core
core_silicate = q.load_material("data/eps_Sil_Oss92.dat", ["wlen", "real_m", "im_m"])
core_silicate.plot_ref_index("ref_index_core.png")

# load and plot mantle
mantle = q.load_material("data/eps_H93.dat", ["wlen", "real_m", "im_m"], units="1/cm")
mantle.plot_ref_index("ref_index_mantle.png")

# compute and plot core opacity
# core.load_kappa("data/kappa_Oss94_naked.dat")
# core.add_plot_kappa("core_kappa.png", postfix=" (Ossenkopf+1994)")


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
#core_silicate.add_plot_kappa("kappa.png", postfix=" (Si Oss93)")


# load amorphous carbon, compute opacity. and plot
core_carbon = q.load_material("data/eps_carb_P93.dat", ["wlen", "real_m", "im_m"])
core_carbon.dust.rho_bulk = 2e0
core_carbon.compute_kappa()
#core_carbon.add_plot_kappa("kappa.png", postfix=" (aC P93)")
# frac = np.array([1e1**(7.55-12.)*28., 1e1**(8.56-12.)*12.*0.2])
frac = np.array([3.63e-5*28., 3.3e-4*12.*0.6])

# frac = np.array([1e0, 0.6])

norm = core_carbon.dust.mass_normalization * frac[1] \
       + core_silicate.dust.mass_normalization * frac[0]

merged = q.merge_kappa([core_silicate, core_carbon], frac, mass_normalization=norm)
merged.add_plot_kappa("kappa.png", postfix=" (merged, bare)",
                      linestyle="-")

# make composite materials (Si+ice and aC+ice)
composite_silicate = q.make_optical([core_silicate, mantle])
composite_carbon = q.make_optical([core_carbon, mantle])

# compute kappa for different shell thicknesses
for vratio in [0.5, 4.5]:
    aratio = ((vratio + 1e0)**(1./3.) - 1e0)
    composite_silicate.dust.aratio = aratio
    composite_silicate.compute_kappa()
    #composite_silicate.add_plot_kappa("kappa.png", postfix=" (Si+ice, ratio %.1f)" % vratio)

    composite_carbon.dust.aratio = aratio
    composite_carbon.dust.rho_bulk = 2e0
    composite_carbon.compute_kappa()
    #composite_carbon.add_plot_kappa("kappa.png", postfix=" (aC+ice, ratio %.1f)" % vratio)

    merged = q.merge_kappa([composite_silicate, composite_carbon], frac, mass_normalization=norm)
    merged.add_plot_kappa("kappa.png", postfix=" (merged, ratio %.1f)" % vratio,
                          linestyle="-")

