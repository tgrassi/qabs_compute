from qabsmanager import QabsManager

# create object
q = QabsManager()


# load core data from file
core = q.load_material("data/eps_carb_P93.dat", labels=["wlen", "real_m", "im_m"])

# compute opacity for core material
core.compute_kappa()

# print some info to screen
q.report()

# save opacity to plot
core.plot_kappa("kappa_01.png")
