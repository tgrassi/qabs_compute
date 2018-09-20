import numpy as np


class Qabs_utils:
    import matplotlib.pyplot as plt

    # ******************
    def __init__(self):
        self.data = dict()
        self.data_ref = dict()
        self.data_q = dict()
        self.ptype_default = None
        self.data_k = dict()

        self.constants = {"clight": 2.99792458e10,
                          "hplanck_eV*s": 4.135667662e-15}

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
    def load_eps(self, fname, labs=None, reverse=True):
        import sys

        # empty data
        self.data = dict()

        # default labs input
        default_labs = ["wlen", "real_eps1", "im_eps", "real_m1", "im_m"]

        # use default is labs is not set
        if labs is None:
            labs = default_labs

        # init data dictionary
        data = {lab: [] for lab in labs}
        # loop on file to read
        for row in open(fname, "rb"):
            srow = row.strip().replace("\t", " ")
            if row.startswith("#") or srow == "":
                continue
            arow = [float(x) for x in srow.split(" ") if x != ""]
            if len(labs) != len(arow):
                sys.exit("ERROR: file " + fname +
                         " has more columns than the list provided in labs argument!")
            for ii, lab in enumerate(labs):
                data[lab].append(arow[ii])

        # reverse data if needed
        if reverse:
            irev = -1
        else:
            irev = 1
        data = {lab: np.array(x[::irev]) for lab, x in data.iteritems()}

        # convert to Hz
        data["freq"] = self.constants["clight"] / (data["wlen"] * 1e-4)

        # convert to eV
        data["energy_eV"] = self.constants["hplanck_eV*s"] * data["freq"]

        # find real_m from real_m1 if needed
        if "real_m1" in data:
            print "NOTE: real_m computed from real_m1"
            data["real_m"] = data["real_m1"] + 1e0

        # find real_m from real_eps1 if needed
        if "real_eps1" in data:
            print "NOTE: real_eps computed from real_eps1"
            data["real_eps"] = data["real_eps1"] + 1e0

        # compute refractive index from eps
        e1 = data["real_eps"]
        e2 = data["im_eps"]
        data["im_m_computed"] = ((-e1 + np.sqrt(e1**2 + e2**2)) / 2.)**0.5
        data["real_m_computed"] = e2 / 2. / data["im_m_computed"]

        # compute im_m from eps if missing
        if "im_m" not in data:
            print "NOTE: im_m computed from eps"
            data["im_m"] = data["im_m_computed"]

        # compute real_m from eps if missing
        if "real_m" not in data:
            print "NOTE: real_m computed from eps"
            data["real_m"] = data["real_m_computed"]

        # compute im_m1 from im_m if missing
        if "real_m1" not in data:
            print "NOTE: real_m1 computed from real_m"
            data["real_m1"] = data["real_m"] - 1e0

        # remove dummy if present
        if "dummy" in data:
            data["dummy"] = None

        print "Loading from " + fname + " done!"

        # copy data to class attribute
        self.data = data

    # ************************
    # perform a benchmark
    def benchmark(self):
        # Qabs benchmark
        self.load_eps("eps_Sil.dat", labs=["wlen", "real_eps1", "im_eps", "real_m1", "im_m"])
        self.load_qref("Sil_21_1e3.dat")
        self.compute_q(1e-7)
        self.plot(what=["qabs", "qabs_ref", "qsca", "qsca_ref"], fname="benchmark_q.png",
                  styles=["-", ":", "-", ":"])

        # computed refractive index
        self.plot(what=["real_m", "im_m", "real_m_computed", "im_m_computed"],
                  fname="benchmark_m.png", styles=["-", "-", ":", ":"])

        # compute opacity
        self.compute_kappa()
        self.plot(what=["kappa"], fname="benchmark_kappa.png")

    # *****************
    # compute opacity (cross section per dust mass), cm2/g
    def compute_kappa(self, amin=5e-7, amax=2.5e-5, pexp=-3.5, rho_bulk=2.9):

        pexp4 = pexp + 4.
        ngrid = 30  # number of size grid points

        # size range, cm
        arange = np.logspace(np.log10(amin), np.log10(amax), ngrid)

        qdata = []
        # loop on sizes
        for ii, asize in enumerate(arange):
            print round(ii*1e2 / (ngrid - 1), 1), "%"
            self.compute_q(asize)
            qabs = self.data_q["qabs"]
            qdata.append(qabs)
        qdata = np.array(qdata).T

        # size distribution
        phi = arange**pexp

        # normalized mass
        cnorm = 4./3. * np.pi * rho_bulk * (amax**pexp4 - amin**pexp4) / pexp4
        icnorm = 1e0 / cnorm

        kappa = []
        # loop to compute integral
        for q in qdata:
            kappa.append(np.pi * np.trapz(q * arange**2*phi, arange) * icnorm)

        # copy to attribute
        self.data_k["wlen"] = self.data_q["wlen"]
        self.data_k["freq"] = self.data_q["freq"]
        self.data_k["kappa"] = kappa

    # ************************
    # load qabs data for comparison
    def load_qref(self, fname):

        # empty data
        self.data_ref = dict()

        labs = ["wlen", "qabs_ref", "qsca_ref", "gcos_ref"]
        data_qref = {lab: [] for lab in labs}
        for row in open(fname, "rb"):
            srow = row.strip()
            if row.startswith("#") or srow == "":
                continue
            arow = [float(x) for x in srow.split(" ") if x != ""]
            for ii, lab in enumerate(labs):
                data_qref[lab].append(arow[ii])
        self.data_ref = {lab: np.array(x) for lab, x in data_qref.iteritems()}

    # ************************
    # compute qabs and qsca
    # data is the data structure loaded by load_eps
    def compute_q(self, asize, nang=100):
        from bhmie import bhmie

        # local copy of data to work with
        data = self.data

        data_q = {"qabs": [], "qsca": [], "qbak": []}
        # loop on wavelengths to call bhmie
        for ii, wlen in enumerate(data["wlen"]):
            # non-dimensional parameter where 1e4 is to convert asize to micron
            xr = 2. * np.pi / wlen * asize * 1e4
            # imaginary refractive index as complex number
            refrel = complex(data["real_m"][ii], data["im_m"][ii])
            # call bhmie that does the magic
            s1, s2, qext, qsca, qbak, gsca = bhmie(xr, refrel, nang)
            # store data
            data_q["qabs"].append(qext-qsca)
            data_q["qsca"].append(qsca)
            data_q["qbak"].append(qbak)
        # stores wavelengths
        data_q["wlen"] = data["wlen"]
        # convert to Hz
        data_q["freq"] = self.constants["clight"] / (data_q["wlen"] * 1e-4)

        # convert to numy arrays
        self.data_q = {lab: np.array(x) for lab, x in data_q.iteritems()}

    # ******************************
    def plot(self, xref="wlen", what=None, fname="output.png", ptype=plt.loglog, styles=None):
        import matplotlib.pyplot as plt
        import sys

        # set default if empty
        if what is None:
            what = ["qsca", "qabs"]

        if styles is None:
            styles = ["-"] * len(what)

        # use plot type default is set
        if self.ptype_default is not None:
            ptype = self.ptype_default

        # if ptype is a string eval to matplotlib function
        if type(ptype) is str:
            if ptype in ["linear", "linlin"]:
                ptype = "plot"
            ptype = eval("plt." + ptype)

        plt.clf()
        # check where to take data from
        for ii, wh in enumerate(what):
            if wh in self.data.keys():
                data_plot = self.data
            elif wh in self.data_q.keys():
                data_plot = self.data_q
            elif wh in self.data_ref.keys():
                data_plot = self.data_ref
            elif wh in self.data_k.keys():
                data_plot = self.data_k
            else:
                sys.exit("ERROR: unknonw lab " + wh + " to plot")

            ptype(data_plot[xref], data_plot[wh], label=wh, linestyle=styles[ii])

        if xref == "wlen":
            xlabel = "$\\lambda / \\mu$m"
        elif xref == "freq":
            xlabel = "$\\nu$ / Hz"
        elif xref == "energy_eV":
            xlabel = "$E$ / eV"
        else:
            sys.exit("ERROR: unknonw xref " + xref + " for the xlabel")

        # add some ornaments to plot
        plt.xlabel(xlabel)
        plt.ylabel(" or ".join(what))
        plt.legend(loc="best")
        plt.savefig(fname)

        print "Plot saved to " + fname
