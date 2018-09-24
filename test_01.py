from utils import QabsManager

# create object
q = QabsManager()

# load data from file with given format
core = q.load_material("eps_Sil.dat", name="core")

core.plot_ref_index("ref_index.png")

# load data for mantle/coating material
mantle = q.load_material("water_ice.dat", labs=["wlen", "real_m", "im_m"], name="mantle")

core.compute_kappa()
core.save_kappa("test.dat")


asdjlksd

# write a report with some information on the loaded data
q.report()

# create new composite material
icy = q.make_optical([core, mantle])

core.compute_q(1e-6)
mantle.compute_q(1e-6)
icy.compute_q([1e-6, 1.2e-6])

icy.plot_q()

icy.load_kappa("test.dat")

# icy.dust.alayer = 1e-7
# icy.compute_kappa()
# icy.save_kappa("test.dat")

icy.plot_kappa()

dadasdas

icy.dust.alayer = 1e-6
icy.compute_kappa()
icy.dust.add_plot_kappa("output.png")
