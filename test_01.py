from utils import Qabs_utils

# create object
q = Qabs_utils()

# load data from file with given format
q.load_eps("eps_Sil.dat")

# compute qabs and qsca for give grain size
asize = 1e-7  # cm
q.compute_q(asize)

# plot eps
q.plot(what=["real_eps", "im_eps"], fname="test_01_eps.png")

# plot refractive index
q.plot(what=["real_m", "im_m"], fname="test_01_m.png")

# plot qabs and qsca
q.plot(what=["qabs", "qsca"], fname="test_01_q.png")

