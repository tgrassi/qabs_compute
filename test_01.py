from utils import QabsManager

# create object
q = QabsManager()

# load core data from file
core = q.load_material("data/eps_Sil.dat")

# plot refractive index
core.plot_ref_index("ref_index_core.png")

# compute opacity for core material
core.compute_kappa()

# store core opacity to file
core.save_kappa("kappa_core.out")


# load data for mantle/coating material
mantle = q.load_material("data/water_ice.dat", labs=["wlen", "real_m", "im_m"])

# create new composite material, core and mantle
icy = q.make_optical([core, mantle])


# write a report with some information on the loaded data
q.report()


# compute Qabs for a given size
core.compute_q(1e-6)
mantle.compute_q(1e-6)

# plot Qabs
core.plot_q("qabs_core.png")
mantle.plot_q("qabs_mantle.png")


# compute Qabs for the core+mantle system
icy.compute_q([1e-6, 1.2e-6])

# plot Qabs
icy.plot_q("qabs_icy.png")

# set layer thickness
icy.dust.alayer = 1e-7

# print some info on the dust used
icy.dust.report()

# compute opacity for the compound material
icy.compute_kappa()

# store to file
icy.save_kappa("kappa_icy.out")

# plot opacity to file
icy.plot_kappa("kappa_icy_1e7.png")

# change layer thickness
icy.dust.alayer = 1e-6

# compute opacity
icy.compute_kappa()

# plot opacity to file
icy.dust.add_plot_kappa("kappa_icy_1e6.png")
