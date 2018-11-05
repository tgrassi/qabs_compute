

class Material:

    # ******************
    def __init__(self, fname=None, labs=None, name=None, units="micron"):

        self.data = dict()

        if fname is not None:
            self.name = name
            self.load(fname, labs, units)

    # ***********************
    # add a single point value
    def add_value(self, wlen, ref_index):
        import numpy as np

        if "wlen" not in self.data:
            self.data["wlen"] = []
            self.data["real_m"] = []
            self.data["im_m"] = []

        self.data["wlen"].append(wlen)
        self.data["real_m"].append(np.real(ref_index))
        self.data["im_m"].append(np.imag(ref_index))

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
    # Returns a dictionary of numpy arrays with keys as labs
    def load(self, fname, labs, units):
        import sys
        import numpy as np
        from utility import hz_to_ev, micron_to_hz, ev_to_micron, wavenumber_to_micron

        print "Loading from " + fname + "..."

        # wavelength is mandatory
        if "wlen" not in labs:
            sys.exit("ERROR: wlen label is needed when loading data!")

        allowed_labs = ["wlen", "real_m", "im_m", "real_m1", "real_eps",
                        "im_eps", "real_eps1"]

        for lab in labs:
            if lab not in allowed_labs:
                print "ERROR: found unknown label " + lab
                print " allowed labels are", allowed_labs
                sys.exit()

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

        # convert to units if necessary
        if units.lower() == "micron":
            print "Found units " + units
        elif units.lower() == "ev":
            print "Found units " + units + ", converting to micron"
            data["wlen"] = ev_to_micron(data["wlen"])
        elif units.lower() == "1/cm":
            print "Found units " + units + ", converting to micron"
            data["wlen"] = wavenumber_to_micron(data["wlen"])
        elif units.lower() == "cm":
            print "Found units " + units + ", converting to micron"
            data["wlen"] = np.array(data["wlen"]) * 1e4
        else:
            sys.exit("ERROR: unknown units " + units)

        # reverse data if needed
        if data["wlen"][0] > data["wlen"][1]:
            print "reversing data..."
            irev = -1
        else:
            irev = 1

        # lists to numpy arrays
        data = {lab: np.array(x[::irev]) for lab, x in data.iteritems()}

        # convert to Hz
        data["freq"] = micron_to_hz(data["wlen"])

        # convert to eV
        data["energy_eV"] = hz_to_ev(data["freq"])

        # find real_m from real_m1 if needed
        if "real_m1" in data and "real_m" not in data:
            print "NOTE: real_m computed from real_m1"
            data["real_m"] = data["real_m1"] + 1e0

        # find real_eps from real_eps1 if needed
        if "real_eps1" in data and "real_eps" not in data:
            print "NOTE: real_eps computed from real_eps1"
            data["real_eps"] = data["real_eps1"] + 1e0

        # compute real_m and im_m from eps
        if "real_eps" in data and "im_eps" in data:
            # compute refractive index from eps
            e1 = data["real_eps"]
            e2 = data["im_eps"]
            data["im_m_computed"] = ((-e1 + np.sqrt(e1 ** 2 + e2 ** 2)) / 2.) ** 0.5
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

        # compute eps if missing from im_m and real_m
        if "real_eps" not in data:
            data["real_eps"] = data["real_m"]**2 - data["im_m"]**2
            data["im_eps"] = 2e0 * data["real_m"] * data["im_m"]

        # remove dummy if present
        if "dummy" in data:
            data["dummy"] = None

        # store file name
        data["fname"] = fname

        print "Loading from " + fname + " done!"

        # copy data to attribute
        self.data = data

    # ***********************
    def extrapolate(self, wmax, nlast=30):
        import numpy as np
        import sys
        from scipy.optimize import curve_fit
        from lmfit import minimize, Parameters, Parameter

        def ff(x, fl, wl2, gamma):
            clight = 3e10
            echarge = 4.80320425e-10
            me = 9.10938e-28
            pre = 2e0 * echarge**2 / me / clight
            return fl * wl2 * x**2 / (2e0 * np.pi * clight * (x**2 - wl2)
                                      + 1j * gamma * wl2 * x)

        def residuals(paramz, x, data):
            fl = paramz['fl'].value
            wl2 = paramz['wl2'].value
            gamma = paramz['gamma'].value

            diff = ff(x, fl, wl2, gamma) - data
            # The only change required is to use view instead of abs.
            return diff.view(np.float)

        params = Parameters()
        params.add("fl", value=1e15, min=1e9, max=1e20)
        params.add("wl2", value=1e-4, min=1e-10, max=1e-1)
        params.add("gamma", value=1e14, min=1e8, max=1e20)

        th = self.data["wlen"][-nlast:] * 1e-4
        ydata = (self.data["real_m"][-nlast:] + 1j * self.data["im_m"][-nlast:])**2 - 1e0

        result = minimize(residuals, params, args=(th, ydata), method="leastsq")
        fl = result.params["fl"].value
        wl2 = result.params["wl2"].value
        gamma = result.params["gamma"].value

        print result.params

        zdata = ff(th, fl, wl2, gamma)

        import matplotlib.pyplot as plt
        plt.clf()
        plt.plot(th, ydata)
        plt.plot(th, zdata)
        plt.show()

        askljdlkas

        def fe1(x2, wp2, wo2, gamma):
            return 1e0 + wp2 * (wo2 - x2) / ((wo2 - x2)**2 + gamma**2 * x2)

        def fe2(x, wp2, wo2, gamma):
            return wp2 * gamma * x / ((wo2 - x**2)**2 + gamma**2 * x**2)

        step = 1
        c1 = c2 = None
        s1 = s2 = None
        for _ in range(10):
            c2, s2 = curve_fit(fe2, 1e0 / self.data["wlen"][-nlast::step],
                               self.data["im_eps"][-nlast::step], p0=c1)
            c1, s1 = curve_fit(fe1, 1e0 / self.data["wlen"][-nlast::step]**2,
                               self.data["real_eps"][-nlast::step], p0=c2)
        if np.inf in s1 or np.inf in s2:
            sys.exit("ERROR: infinity in covariance matrix!")

        wmax_data = self.data["wlen"][-1]
        for wlen in np.linspace(wmax_data, wmax, 20):
            e1 = fe1(1e0 / wlen**2, c1[0], c1[1], c1[2])
            e2 = fe2(1e0 / wlen, c1[0], c1[1], c1[2])
            m2 = ((-e1 + np.sqrt(e1 ** 2 + e2 ** 2)) / 2.) ** 0.5
            m1 = e2 / 2. / m2

            self.data["wlen"] = np.append(self.data["wlen"], wlen)
            self.data["real_m"] = np.append(self.data["real_m"], m1)
            self.data["im_m"] = np.append(self.data["im_m"], m2)
            self.data["real_eps"] = np.append(self.data["real_eps"], e1)
            self.data["im_eps"] = np.append(self.data["im_eps"], e2)

    # ******************
    # save opacity data to file
    def save_refractive_index(self, fname):

        print "Saving material to file " + fname
        fout = open(fname, "w")
        fout.write("# wlen(micron), Re(eps), Im(eps), Re(m), Im(m)\n")
        for ii, wlen in enumerate(self.data["wlen"]):
            row = [str(self.data[k][ii]) for k in ["wlen", "real_eps", "im_eps", "real_m", "im_m"]]
            fout.write(" ".join(row) + "\n")
        fout.close()
