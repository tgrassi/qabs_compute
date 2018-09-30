from utils import QabsManager
import numpy as np

# create object
q = QabsManager()

# load and plot core
core = q.load_material("data/eps_Sil_Oss92.dat", ["wlen", "real_m", "im_m"])
core.plot_ref_index("ref_index_core.png")

# load and plot mantle
mantle = q.load_material("data/eps_H93.dat", ["wlen", "real_m", "im_m"], units="1/cm")
mantle.plot_ref_index("ref_index_mantle.png")

# compute and plot core opacity
# core.load_kappa("data/kappa_Oss94_naked.dat")
# core.add_plot_kappa("core_kappa.png", postfix=" (Ossenkopf+1994)")


# compute core opacity and plot
q.clear_plots()
core.compute_kappa()
core.add_plot_kappa("kappa.png", postfix=" (Si Oss93)")


# load amorphous carbon, compute opacity. and plot
core_carbon = q.load_material("data/eps_carb_P93.dat", ["wlen", "real_m", "im_m"])
core_carbon.compute_kappa()
core_carbon.add_plot_kappa("kappa.png", postfix=" (aC P93)")


# make composite materials (Si+ice and aC+ice)
composite_silicate = q.make_optical([core, mantle])
composite_carbon = q.make_optical([core_carbon, mantle])

for vratio in [0.5, 4.5]:
    aratio = ((vratio + 1e0)**(1./3.) - 1e0)
    # compute kappa for different shell thicknesses
    composite_silicate.dust.aratio = aratio
    composite_silicate.compute_kappa()
    composite_silicate.add_plot_kappa("kappa.png", postfix=" (Si+ice, ratio %.1f)" % vratio)

    composite_carbon.dust.aratio = aratio
    composite_carbon.compute_kappa()
    composite_carbon.add_plot_kappa("kappa.png", postfix=" (aC+ice, ratio %.1f)" % vratio)


