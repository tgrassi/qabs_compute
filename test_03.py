from utils import Qabs_utils

# create object
q = Qabs_utils()

# load data from file with given format
q.load_eps("eps_CO.dat", labs=["wlen", "real_eps1", "im_eps"])

# compute qabs and qsca for give grain size
asize = 1e-7  # cm
q.compute_q(asize)

# plot eps
q.plot(xref="freq", what=["real_eps1"], fname="test_03_reps.png", ptype="plot")
q.plot(xref="freq", what=["im_eps"], fname="test_03_ieps.png", ptype="plot")


# plot refractive index
q.plot(xref="freq", what=["real_m"], fname="test_03_rm.png", ptype="plot")
q.plot(xref="freq", what=["im_m"], fname="test_03_im.png", ptype="plot")


# plot qabs and qsca
q.plot(xref="freq", what=["qabs", "qsca"], fname="test_03_q.png", ptype="plot")














