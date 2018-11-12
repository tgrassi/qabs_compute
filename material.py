

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
            print "NOTE: real_m_computed computed from real_eps and im_eps"
            print "NOTE: im_m_computed computed from real_eps and im_eps"
            # compute refractive index from eps
            e1 = data["real_eps"]
            e2 = data["im_eps"]
            data["im_m_computed"] = ((-e1 + np.sqrt(e1 ** 2 + e2 ** 2)) / 2.) ** 0.5
            data["real_m_computed"] = e2 / 2. / data["im_m_computed"]

        # compute im_m from eps if missing
        if "im_m" not in data:
            print "NOTE: im_m computed from im_eps_computed"
            data["im_m"] = data["im_m_computed"]

        # compute real_m from eps if missing
        if "real_m" not in data:
            print "NOTE: real_m computed from real_eps_computed"
            data["real_m"] = data["real_m_computed"]

        # compute im_m1 from im_m if missing
        if "real_m1" not in data:
            print "NOTE: real_m1 computed from real_m"
            data["real_m1"] = data["real_m"] - 1e0

        # compute eps if missing from im_m and real_m
        if "real_eps" not in data:
            print "NOTE: real_eps computed from real_m and im_m"
            print "NOTE: im_eps computed from real_m and im_m"
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
    def extrapolate(self, wmax, nlast=70):
        import numpy as np
        import sys
        from scipy.optimize import curve_fit
        from lmfit import minimize, Parameters, Parameter

        def ffw(x, c_wp, c_wo, c_gamma, c_off):
            fun = c_off + c_wp**2 / (c_wo**2 - x**2 - 1j * x * c_gamma)
            return fun

        def residualsi(paramz, x, data):
            c_wp = paramz['wp'].value
            c_wo = paramz['wo'].value
            c_gamma = paramz['gamma'].value
            c_off = paramz['off'].value

            diff = ffw(x, c_wp, c_wo, c_gamma, c_off) - data

            return diff.view(np.float32)

        clight = 3e10
        th = self.data["wlen"][-nlast:] * 1e-4
        th = 2e0 * np.pi * clight / th
        #ydata = (self.data["real_m"][-nlast:] + 1j * self.data["im_m"][-nlast:])**2 - 1e0
        ydata = self.data["real_eps"][-nlast:] + 1j * self.data["im_eps"][-nlast:]

        params = Parameters()
        params.add("wp", value=1e1**13.47, min=1e1**13.35, max=1e1**13.55)
        params.add("wo", value=1e1**13.58, min=1e1**13.5, max=1e1**13.7)
        params.add("gamma", value=1e1**13.32, min=1e1**13.25, max=1e1**13.4)
        params.add("off", value=2.48, min=2.2, max=2.6)

        result = minimize(residualsi, params, args=(th, ydata), method="leastsq", nan_policy='omit')
        wp = 1e1**13.35 #result.params["wp"].value
        wo = 1e1**13.51 # result.params["wo"].value
        gamma = 1e1**12.92 #result.params["gamma"].value
        off = 2.52 #result.params["off"].value

        for k in result.params.keys():
            print result.params[k]
            print np.log10(result.params[k].value)

        zdata = ffw(th, wp, wo, gamma, off)

        import matplotlib.pyplot as plt
        plt.clf()
        plt.plot(th, ydata.real)
        plt.plot(th, zdata.real)
        plt.show()

        plt.clf()
        plt.plot(th, ydata.imag)
        plt.plot(th, zdata.imag)
        plt.show()

        wmax_data = self.data["wlen"][-1]
        for wlen in np.linspace(wmax_data, wmax, 20):
            ee = ffw(wlen * 1e-4, wp, wo, gamma, off)
            e1 = np.real(ee)
            e2 = 0e0

            cconj = np.sqrt(e1 ** 2 + e2 ** 2)
            m1 = np.sqrt((cconj + e1) / 2e0)
            m2 = np.sqrt((cconj - e1) / 2e0)

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

    # ***********************
    # add an spherical impurities with given volume filling factor.
    # optical properties of the impurities are taken from material_impurity
    def add_impurity(self, material_impurity, volume_filling_factor):
        from scipy.interpolate import interp1d
        from numpy import sqrt

        # interpolate real and imaginary dielectric functions
        feps1 = interp1d(material_impurity.data["wlen"], material_impurity.data["real_eps"])
        feps2 = interp1d(material_impurity.data["wlen"], material_impurity.data["im_eps"])

        # alias name for the volume filling factor
        vff = volume_filling_factor
        for ii, wlen in enumerate(self.data["wlen"]):
            # complex impurities material dielectric
            eps = complex(feps1(wlen), feps2(wlen))
            # complex matrix material dielectric
            eps_m = complex(self.data["real_eps"][ii], self.data["im_eps"][ii])

            # Bruggemann rule, see Bohren+Huffmann 1983, eq. 8.50
            g = vff * (eps - eps_m) / (eps + 2e0 * eps_m)
            eps_av = eps_m * (1e0 + 3e0 * g / (1e0 - g))

            # set new dielectric
            self.data["real_eps"][ii] = eps1 = eps_av.real
            self.data["im_eps"][ii] = eps2 = eps_av.imag
            # compute refractive index from dielectric
            cconj = sqrt(eps1**2 + eps2**2)
            self.data["real_m"][ii] = sqrt((cconj + eps1) / 2e0)
            self.data["im_m"][ii] = sqrt((cconj - eps1) / 2e0)
