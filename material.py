

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
    #  wlen: wavelength in micron
    #  real_eps1: real part of epsilon-1
    #  real_eps: real part of epsilon
    #  im_eps: imaginary part of epsilon
    #  real_m1: real part of refractive index-1
    #  real_m: real part of refractive index
    #  im_m: imaginary part of refractive index
    # Note: you only need wlen, real_m (or real_m1), and im_m
    # Returns a dictionary of numpy arrays with keys as labs
    def load(self, fname, labs, units):
        import sys
        import numpy as np
        from utility import hz_to_ev, micron_to_hz, ev_to_micron, wavenumber_to_micron, thz_to_micron

        print("Loading from " + fname + "...")

        # wavelength is mandatory
        if "wlen" not in labs:
            sys.exit("ERROR: wlen label is needed when loading data!")

        allowed_labs = ["wlen", "real_m", "im_m", "real_m1", "real_eps",
                        "im_eps", "real_eps1", "alpha", "dummy"]

        for lab in labs:
            if lab not in allowed_labs:
                print("ERROR: found unknown label " + lab)
                print(" allowed labels are", allowed_labs)
                sys.exit()

        # init data dictionary
        data = {lab: [] for lab in labs}
        # loop on file to read
        for row in open(fname):
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
            print("Found units " + units)
        elif units.lower() == "ev":
            print("Found units " + units + ", converting to micron")
            data["wlen"] = ev_to_micron(data["wlen"])
        elif units.lower() == "1/cm":
            print("Found units " + units + ", converting to micron")
            data["wlen"] = wavenumber_to_micron(data["wlen"])
        elif units.lower() == "cm":
            print("Found units " + units + ", converting to micron")
            data["wlen"] = np.array(data["wlen"]) * 1e4
        elif units.lower() == "thz":
            print("Found units " + units + ", converting to micron")
            data["wlen"] = thz_to_micron(data["wlen"])
        else:
            sys.exit("ERROR: unknown units " + units)

        # reverse data if needed
        if data["wlen"][0] > data["wlen"][1]:
            print("reversing data...")
            irev = -1
        else:
            irev = 1

        # lists to numpy arrays
        data = {lab: np.array(x[::irev]) for lab, x in data.items()}

        # convert to Hz
        data["freq"] = micron_to_hz(data["wlen"])

        # convert to eV
        data["energy_eV"] = hz_to_ev(data["freq"])

        # find real_m from real_m1 if needed
        if "real_m1" in data and "real_m" not in data:
            print("NOTE: real_m computed from real_m1")
            data["real_m"] = data["real_m1"] + 1e0

        # find real_eps from real_eps1 if needed
        if "real_eps1" in data and "real_eps" not in data:
            print("NOTE: real_eps computed from real_eps1")
            data["real_eps"] = data["real_eps1"] + 1e0

        # compute real_m and im_m from eps
        if "real_eps" in data and "im_eps" in data:
            print("NOTE: real_m_computed computed from real_eps and im_eps")
            print("NOTE: im_m_computed computed from real_eps and im_eps")
            # compute refractive index from eps
            e1 = data["real_eps"]
            e2 = data["im_eps"]
            data["im_m_computed"] = ((-e1 + np.sqrt(e1 ** 2 + e2 ** 2)) / 2.) ** 0.5
            data["real_m_computed"] = e2 / 2. / data["im_m_computed"]

        if "alpha" in data and "im_eps" not in data and "real_m" in data:
            print("NOTE: im_eps computed from alpha")
            data["im_m"] = data["alpha"] * data["wlen"] / 1e4 / 2e0 / np.pi
            data["real_eps"] = data["real_m"]**2 - data["im_m"]**2
            data["im_eps"] = 2e0 * data["real_m"] * data["im_m"]
            data["im_m_computed"] = data["im_m"]

        # compute im_m from eps if missing
        if "im_m" not in data:
            print("NOTE: im_m computed from im_eps_computed")
            data["im_m"] = data["im_m_computed"]

        # compute real_m from eps if missing
        if "real_m" not in data:
            print("NOTE: real_m computed from real_eps_computed")
            data["real_m"] = data["real_m_computed"]

        # compute im_m1 from im_m if missing
        if "real_m1" not in data:
            print("NOTE: real_m1 computed from real_m")
            data["real_m1"] = data["real_m"] - 1e0

        # compute eps if missing from im_m and real_m
        if "real_eps" not in data:
            print("NOTE: real_eps computed from real_m and im_m")
            print("NOTE: im_eps computed from real_m and im_m")
            data["real_eps"] = data["real_m"]**2 - data["im_m"]**2
            data["im_eps"] = 2e0 * data["real_m"] * data["im_m"]

        # remove dummy if present
        if "dummy" in data:
            data["dummy"] = None

        # store file name
        data["fname"] = fname


        print("Loading from " + fname + " done!")

        # copy data to attribute
        self.data = data

    # ***********************
    # this function extrapolates to wmax (micron) the imaginary part of the eps
    # assuming a 1/wlen profile for the last nlast data point. The real part is
    # computed using Kramers-Kronig on the whole domain.
    # Note, that the extrapolation is pretty arbitrary and might lead to
    # unexpected results. Use with caution.
    def extrapolate(self, wmax, nlast):
        import numpy as np
        from scipy.optimize import curve_fit
        from numpy import log10
        import sys

        clight = 3e10  # cm/s
        th = self.data["wlen"][-nlast:] * 1e-4  # cm
        th = 2e0 * np.pi * clight / th  # Hz

        # eps from data, get only the final nlast data
        ydata = self.data["real_eps"][-nlast:] + 1j * self.data["im_eps"][-nlast:]

        # only include eps with strictly positive imaginary part
        th = [x for ii, x in enumerate(th) if ydata.imag[ii] > 0e0]
        ydata = [x for x in ydata.imag if x > 0e0]

        # function to fit is linear in log space (power-law)
        def fw(xx, a, b):
            return a * xx + b

        # fit in log space
        p, sm = curve_fit(fw, log10(th), log10(ydata))

        # compute fit data to compare with original data
        zdata = 1e1**fw(log10(th), p[0], p[1])

        import matplotlib.pyplot as plt

        plt.clf()
        plt.plot(th, ydata, marker="o", label="data")
        plt.plot(th, zdata, label="fit")
        plt.xlabel("$\omega$ / Hz")
        plt.ylabel("Im(eps)")
        plt.legend(loc="best")
        png_file_fit = "extrapolation_fit_img_eps.png"
        plt.savefig(png_file_fit)
        print("Imaginary part fit saved to " + png_file_fit)

        # number of point to extrapolate after last data point
        npoints = 200

        from copy import copy
        # copy data for later comparison
        wlen_old = copy(self.data["wlen"])
        real_eps_old = copy(self.data["real_eps"])
        im_eps_old = copy(self.data["im_eps"])

        # take only the final nlast data
        self.data["wlen"] = self.data["wlen"][:-nlast]
        self.data["real_eps"] = self.data["real_eps"][:-nlast]
        self.data["im_eps"] = self.data["im_eps"][:-nlast]
        self.data["real_m"] = self.data["real_m"][:-nlast]
        self.data["im_m"] = self.data["im_m"][:-nlast]

        # get last data point from where to extrapolate
        wmax_data = self.data["wlen"][-1]

        # sequence of points to extrapolate
        wextrap = np.linspace(wmax_data, wmax, npoints)

        # extrapolate imaginary part
        e2_extrap = []
        for wlen in wextrap:
            e2 = 1e1 ** fw(log10(2e0 * np.pi * clight / (wlen * 1e-4)), p[0], p[1])
            e2_extrap.append(e2)

        # compute frequency
        ths = np.append(self.data["wlen"], wextrap) * 1e-4
        ths = 2e0 * np.pi * clight / ths

        # apply Kronig-Kramer
        e1_extrap = []
        for wlen in wextrap:
            # compute frequncy for the integral
            thx = 2e0 * np.pi * clight / (wlen * 1e-4)
            # skip poles
            fint = np.array([xx for ii, xx in enumerate(np.append(self.data["im_eps"], e2_extrap)) if ths[ii] != thx])
            xint = np.array([xx for ii, xx in enumerate(ths) if ths[ii] != thx])

            # reverse array if xdata are descending
            if xint[0] > xint[1]:
                fint = fint[::-1]
                xint = xint[::-1]

            # do KK
            kk = 1e0 + 2e0 / np.pi * np.trapz(fint * xint / (xint ** 2 - thx ** 2), xint)

            # check KK sign
            if kk < 0e0:
                sys.exit("ERROR: negative real part of eps from Kramers-Kronig!")

            # store eps real part
            e1_extrap.append(kk)

        # remove first data computed with KK
        e1_extrap = np.array(e1_extrap[1:])
        e2_extrap = e2_extrap[1:]
        wextrap = wextrap[1:]

        # normalize extrapolation to the last data point
        # e1_extrap += self.data["real_eps"][-1] - e1_extrap[0]

        marker = None

        plt.clf()
        plt.semilogx(wlen_old, real_eps_old, marker=marker, label="data")
        plt.semilogx(wextrap, e1_extrap, marker=marker, label="extrapolation")
        plt.xlabel("$\\lambda$ / $\\mu$m")
        plt.ylabel("Re(eps)")
        png_file_fit = "extrapolation_real_eps.png"
        plt.savefig(png_file_fit)
        print("Real part extrapolation saved to " + png_file_fit)

        plt.clf()
        plt.semilogx(wlen_old, im_eps_old, marker=marker, label="data")
        plt.semilogx(wextrap, e2_extrap, marker=marker, label="extrapolation")
        plt.xlabel("$\\lambda$ / $\\mu$m")
        plt.ylabel("Im(eps)")
        png_file_fit = "extrapolation_img_eps.png"
        plt.savefig(png_file_fit)
        print("Imaginary part extrapolation saved to " + png_file_fit)

        # add extrapolate data to material arrays
        for ii, wlen in enumerate(wextrap):
            e1 = e1_extrap[ii]
            e2 = e2_extrap[ii]
            wlen = wextrap[ii]

            self.data["wlen"] = np.append(self.data["wlen"], wlen)
            self.data["real_eps"] = np.append(self.data["real_eps"], e1)
            self.data["im_eps"] = np.append(self.data["im_eps"], e2)

            cconj = np.sqrt(e1 ** 2 + e2 ** 2)
            m1 = np.sqrt((cconj + e1) / 2e0)
            m2 = np.sqrt((cconj - e1) / 2e0)

            self.data["real_m"] = np.append(self.data["real_m"], m1)
            self.data["im_m"] = np.append(self.data["im_m"], m2)

    # ******************
    # save opacity data to file
    def save_refractive_index(self, fname):

        print("Saving material to file " + fname)
        fout = open(fname, "w")
        fout.write("# wlen(micron), Re(eps), Im(eps), Re(m), Im(m)\n")
        for ii, wlen in enumerate(self.data["wlen"]):
            row = ["%e" % self.data[k][ii] for k in ["wlen", "real_eps", "im_eps", "real_m", "im_m"]]
            fout.write(" ".join(row) + "\n")
        fout.close()

    # ***********************
    # add an spherical impurities with given volume filling factor.
    # optical properties of the impurities are taken from material_impurity
    def add_impurity(self, material_impurity_list, volume_filling_factor_list):
        from scipy.interpolate import interp1d
        import numpy as np
        import sys

        # convert argument to list if not
        if type(material_impurity_list) is not list:
            material_impurity_list = [material_impurity_list]

        # convert argument to list if not
        if type(volume_filling_factor_list) is not list:
            volume_filling_factor_list = [volume_filling_factor_list]

        feps1 = []
        feps2 = []
        # interpolate real and imaginary dielectric functions
        for material_impurity in material_impurity_list:
            f1 = interp1d(material_impurity.data["wlen"], material_impurity.data["real_eps"])
            f2 = interp1d(material_impurity.data["wlen"], material_impurity.data["im_eps"])
            feps1.append(f1)
            feps2.append(f2)

        # number of impurities
        n_impurities = len(material_impurity_list)
        if n_impurities > 1:
            sys.exit("ERROR: more than 1 impurity!")
        # np array of volume filling factors
        delta = np.array(volume_filling_factor_list)
        # sum of impurities
        fsum = sum(delta)

        # check amount of impurities
        if fsum > 0.5:
            sys.exit("ERROR: volume filling factor of impurities exceeds 50% of total volume!")

        # loop on frequencies to add impurities
        for ii, wlen in enumerate(self.data["wlen"]):
            print(wlen)
            eps_m = self.data["real_eps"][ii] + 1j * self.data["im_eps"][ii]
            # beta = np.zeros(n_impurities)
            # eps_list = np.zeros(n_impurities, dtype=np.complex)
            # loop on impurities to compute beta
            # for jj in range(n_impurities):
            #    eps_list[jj] = feps1[jj](wlen) + 1j * feps2[jj](wlen)
            #    beta[jj] = 3e0 * eps_m / (eps_list[jj] + 2e0 * eps_m)

            # derived quantities
            # sbete = sum(delta * beta * eps_list)
            # sbet = sum(delta * beta)

            def fbruggeman(eps_a, eps_b, f):
                from numpy import sqrt

                return 0.25 * (-eps_a + 2e0 * eps_b + 3e0 * eps_a * f - 3e0 * eps_b * f +
                               sqrt(8e0 * eps_a * eps_b + (-eps_a + 2e0 * eps_b
                                                           + 3e0 * eps_a * f - 3e0 * eps_b * f)**2))

            eps_impurity = feps1[0](wlen) + 1j * feps2[0](wlen)
            eps_av = fbruggeman(eps_impurity, eps_m, delta[0])

            # average eps
            # eps_av = ((1e0 - fsum) * eps_m + sbete) / (1e0 - fsum + sbet)
            # sfe = sum(delta * eps_list) + (1e0 - sum(delta)) * eps_m
            # sfde = sum(delta / eps_list) + (1e0 - sum(delta)) / eps_m
            # eps_av = 0.5 * (sfe + 1e0 / sfde)

            # set new dielectric
            self.data["real_eps"][ii] = eps1 = eps_av.real
            self.data["im_eps"][ii] = eps2 = eps_av.imag

            # compute refractive index from dielectric
            cconj = np.sqrt(eps1**2 + eps2**2)
            self.data["real_m"][ii] = np.sqrt((cconj + eps1) / 2e0)
            self.data["im_m"][ii] = np.sqrt((cconj - eps1) / 2e0)
