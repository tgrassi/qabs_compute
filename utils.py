import numpy as np

# *******************************
# load refractive index data from a refractive index file named fname
# labs are the column labels, where
# wlen: wavelength in micron
# real_eps1: real part of epsilon-1
# real_eps: real part of epsilon
# im_eps: imaginary part of epsilon
# real_m1: real part of refractive index-1
# real_m: real part of refractive index
# im_m: imaginary part of refractive index
# Note: you only need wlen, real_m (or real_m1), and im_m
# Use reverse if wavelengths are stored as energies (high to low)
# Returns a dictionary of numpy arrays with keys as labs
def load_eps(fname, labs=["wlen", "real_eps1", "im_eps", "real_m1", "im_m"], reverse=False):
    # init data dictionary
    data = {lab:[] for lab in labs}
    # loop on file to read
    for row in open(fname, "rb"):
        srow = row.strip()
        if row.startswith("#") or srow == "":
            continue
        arow = [float(x) for x in srow.split(" ") if x != ""]
        for ii, lab in enumerate(labs):
            data[lab].append(arow[ii])

    # reverse data if needed
    if reverse:
        irev = -1
    else:
        irev = 1
    data = {lab: np.array(x[::irev]) for lab,x in data.iteritems()}

    # find real_m from real_m1 if needed
    if "real_m1" in data:
        data["real_m"] = data["real_m1"] + 1e0

    # find real_m from real_eps1 if needed
    if "real_eps1" in data:
        data["real_eps"] = data["real_eps1"] + 1e0

    # compute im_m from eps if missing
    if "im_m" not in data:
        e1 = data["real_eps"]
        e2 = data["im_eps"]
        data["im_m"] = ((-e1 + np.sqrt(e1**2 + e2**2)) / 2.)**0.5

    # compute real_m from eps if missing
    if "real_m" not in data:
        e2 = data["im_eps"]
        data["real_m"] = e2 / 2. / data["im_m"]

    return data

# ************************
# load qabs data for comparison
def load_qref(fname):
    labs = ["wlen", "qabs", "qsca", "gcos"]
    data_qref = {lab:[] for lab in labs}
    for row in open(fname, "rb"):
        srow = row.strip()
        if row.startswith("#") or srow == "":
            continue
        arow = [float(x) for x in srow.split(" ") if x != ""]
        for ii, lab in enumerate(labs):
            data_qref[lab].append(arow[ii])
    return {lab: np.array(x) for lab,x in data_qref.iteritems()}

# ************************
# compute qabs and qsca
# data is the data structure loaded by load_eps
def get_q(data, asize):
    from bhmie import bhmie
    nang = 1000
    data_q = {"qabs": [], "qsca": [], "qbak":[]}
    # loop on wavelengths to call bhmie
    for ii, wlen in enumerate(data["wlen"]):
        # non-dimensional parameter where 1e4 is to convert asize to micron
        xr = 2. * np.pi / wlen * asize * 1e4
        # imaginary refractive index as complex number
        refrel = complex(data["real_m"][ii], data["im_m"][ii])
        # call bhmie that does the magic
        S1, S2, qext, qsca, qbak, gsca = bhmie(xr, refrel, nang)
        # store data
        data_q["qabs"].append(qext-qsca)
        data_q["qsca"].append(qsca)
        data_q["qbak"].append(qbak)
    # stores wavelengths
    data_q["wlen"] = data["wlen"]
    # convert to numy arrays
    return {lab: np.array(x) for lab,x in data_q.iteritems()}

