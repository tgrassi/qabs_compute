from utils import QabsManager
from optical import Optical
import numpy as np

# create object
q = QabsManager()

opt = Optical(None)

opt.load_kappa("data/kappa_oss94_v00.dat")
opt.add_plot_kappa("kappa.png", postfix=" (ref, ratio 0.0)", linestyle=":")

opt.load_kappa("data/kappa_oss94_v05.dat")
opt.add_plot_kappa("kappa.png", postfix=" (ref, ratio 0.5)", linestyle=":")

opt.load_kappa("data/kappa_oss94_v45.dat")
opt.add_plot_kappa("kappa.png", postfix=" (ref, ratio 4.5)", linestyle=":")

# load core optical properties
core_silicate = q.load_material("data/eps_Sil_Oss92.dat", ["wlen", "real_m", "im_m"])
core_carbon = q.load_material("data/eps_carb_P93.dat", ["wlen", "real_m", "im_m"])

# load mantle optical properties
mantle = q.load_material("eps_CO.dat", labs=["wlen", "real_eps", "im_eps"])

# add mantle to cores
composite_silicate = q.make_optical([core_silicate, mantle])
composite_carbon = q.make_optical([core_carbon, mantle])

core_silicate.compute_kappa()
core_carbon.dust.rho_bulk = 2e0
core_carbon.compute_kappa()

# volume fractions
frac = np.array((1., 0.475862068965517))
frac /= sum(frac)

merged = q.merge_kappa([core_silicate, core_carbon], frac)
merged.add_plot_kappa("kappa.png", postfix=" (merged, ratio 0.0)",
                      linestyle="-")

# compute kappa for different shell thicknesses
for vratio in [0.5, 4.5]:
    # size ratio from volume ratio
    aratio = ((vratio + 1e0)**(1./3.) - 1e0)

    # silicate + mantle opacity
    composite_silicate.dust.aratio = aratio
    composite_silicate.compute_kappa()

    # carbon + mantle opacity
    composite_carbon.dust.aratio = aratio
    composite_carbon.dust.rho_bulk = 2e0  # g/cm3
    composite_carbon.compute_kappa()

    merged = q.merge_kappa([composite_silicate, composite_carbon], frac)
    merged.add_plot_kappa("kappa.png", postfix=" (merged, ratio %.1f)" % vratio,
                          xlim=(1e2, None), ylim=(None, 1e2))











