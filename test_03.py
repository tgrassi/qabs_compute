from qabsmanager import QabsManager

# create object
q = QabsManager()

q.clear_plots()

# load core data from file
core = q.load_material("data/eps_carb_P93.dat", labels=["wlen", "real_m", "im_m"])

# compute opacity for core material
core.compute_kappa()

# plot core material opacity
core.add_plot_kappa("kappa_03.png", postfix="bare")


# load ice data from file
mantle = q.load_material("data/eps_H93.dat", labels=["wlen", "real_m", "im_m"])

# create new material from core and mantle
composite = q.make_optical([core, mantle])

# volume ratio to size ratio
vratio = 0.5
aratio = ((vratio + 1e0) ** (1. / 3.) - 1e0)

# set size ratio silicon/ice material
composite.dust.aratio = aratio

# compute opacity of the composite material
composite.compute_kappa()

# plot opacity
composite.add_plot_kappa("kappa_03.png", postfix="coated", xlim=(5e1, 1e3), ylim=(1, 1e3))
