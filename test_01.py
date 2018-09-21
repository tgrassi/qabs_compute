from utils import Qabs_utils

# create object
q = Qabs_utils()

# load data from file with given format
q.load_eps("eps_Sil.dat")

# load data for mantle/coating material
q.load_eps_coating("water_ice.dat", labs=["wlen", "real_m", "im_m"])

# write a report with the information on the loaded data
q.report()

# compute qabs and qsca for give grain size
asize = 1e-6  # cm
q.compute_q(asize)

# compute qabs and qsca with coating
asize_coat = asize + 3e-7  # add coating, cm
q.compute_q_coating(asize, asize_coat)

# compute opacity
q.compute_kappa()

# compute opacity adding a layer
q.compute_kappa_coating(alayer=1e-6)

# plot eps
q.plot(what=["real_eps", "im_eps"], fname="test_01_eps.png")

# plot refractive index
q.plot(what=["real_m", "im_m"], fname="test_01_m.png", db_types=["standard", "coating"])

# plot qabs and qsca
q.plot(what=["qabs", "qsca"], fname="test_01_q.png", db_types=["standard", "coating"])

# plot opacity
q.plot(what=["kappa"], fname="test_01_k.png", db_types=["standard", "coating"])
