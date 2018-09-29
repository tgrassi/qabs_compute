from utils import QabsManager
import numpy as np

# create object
q = QabsManager()

# load data from file with given format
mantle = q.load_material("eps_CO.dat", labs=["wlen", "real_eps", "im_eps"])

# load core optical properties
core = q.load_material("data/eps_Sil.dat")

# create new material
composite = q.make_optical([core, mantle])


# QABS
# compute qabs and qsca for a given grain size and variable mantle radii
asize_core = 1e-7  # cm
# loop on mantle radius (note: radius not shell thickness)
for asize_mantle in np.linspace(1.1e-7, 3.1e-7, 5):
    print "mantle radius, cm:", asize_mantle
    composite.compute_q([asize_core, asize_mantle])
    # plot
    composite.add_plot_q("composite_q.png", postfix=(" (%.1e cm)" % asize_mantle))


# OPACITY
q.clear_plots()
# compute kappa for different shell thicknesses
for alayer in np.logspace(-7., -5., 3):
    print "tickness, cm:", alayer
    composite.dust.alayer = alayer
    composite.compute_kappa()
    # plot
    composite.add_plot_kappa("composite_kappa.png", postfix=(" (%.1e cm)" % alayer))













