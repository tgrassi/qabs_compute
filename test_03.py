from utils import Qabs_utils

# create object
q = Qabs_utils()

# load data from file with given format
q.load_eps("eps_CO.dat", labs=["wlen", "real_eps1", "im_eps"])

# compute qabs and qsca for give grain size
asize = 1e-7  # cm
q.compute_q(asize)

# set plot default, in this case linear-linear
q.ptype_default = "linlin"

# plot eps
q.plot(xref="energy_eV", what=["real_eps1"], fname="test_03_reps.png")
q.plot(xref="freq", what=["im_eps"], fname="test_03_ieps.png")


# plot refractive index
q.plot(xref="freq", what=["real_m"], fname="test_03_rm.png")
q.plot(xref="freq", what=["im_m"], fname="test_03_im.png")


# plot qabs and qsca
q.plot(xref="freq", what=["qabs", "qsca"], fname="test_03_q.png")














