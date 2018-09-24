

class Dust:

    # ********************
    def __init__(self, optical=None, amin=5e-7, amax=2.5e-5, pexp=-3.5, rho_bulk=2.9,
                 ngrid=30, alayer=None):
        self.optical = optical
        self.amin = amin
        self.amax = amax
        self.pexp = pexp
        self.rho_bulk = rho_bulk
        self.ngrid = ngrid
        self.alayer = alayer
        self.data = dict()

    # ******************
    def plot_kappa(self, fname):
        import matplotlib.pyplot as plt
        plt.clf()
        self.add_plot_kappa(fname)

    # *******************
    def add_plot_kappa(self, fname):
        import matplotlib.pyplot as plt

        print "Plotting kappa to " + fname + "..."

        plt.loglog(self.data["wlen"], self.data["kappa"])
        plt.xlabel("$\\lambda$ / $\\mu$m")
        plt.ylabel("$\\kappa$ / [cm$^2$ g$^{-1}$]")
        plt.savefig(fname)

    # ******************
    def save_kappa(self, fname):

        fout = open(fname, "w")
        for ii, wlen in enumerate(self.data["wlen"]):
            fout.write(str(wlen) + " " + str(self.data["kappa"][ii]) + "\n")
        fout.close()

    # *********************
    def load_kappa(self, fname):

        print "Loading kappa from " + fname + "..."
        self.data = {"wlen": [], "kappa": []}
        for row in open(fname, "rb"):
            srow = row.strip()
            if srow.startswith("#"):
                continue
            if srow == "":
                continue
            wlen, kappa = [float(x) for x in srow.split(" ") if x != ""]
            self.data["wlen"].append(wlen)
            self.data["kappa"].append(kappa)
        print "done"

    # *********************
    def compute_kappa(self, verbose=1):
        import sys

        if len(self.optical.materials) == 1:
            self.compute_kappa_bare(verbose)
        elif len(self.optical.materials) == 2:
            self.compute_kappa_coating(verbose)
        else:
            sys.exit("ERROR: cannot compute kappa with more than 2 materials")

    # *****************
    # compute opacity (cross section per dust mass), cm2/g
    # Assumes MNR dust size distribution with power-law exponent pexp,
    # in the range amin, amax, rho_bulk bulk density. Integral on dust
    # distribution uses ngrid points
    def compute_kappa_bare(self, verbose=1):
        import numpy as np

        # power-law exponent + 4 for integrals
        pexp4 = self.pexp + 4.

        # size range, cm
        arange = np.logspace(np.log10(self.amin), np.log10(self.amax), self.ngrid)

        if verbose > 0:
            print "Computing opacity..."

        qdata = []
        # loop on grain sizes
        for ii, asize in enumerate(arange):
            if verbose > 0:
                print round(ii * 1e2 / (self.ngrid - 1), 1), "%"
            self.optical.compute_q(asize)
            qabs = self.optical.data["qabs"]
            qdata.append(qabs)
        # transpose from (asize, wlen) to (wlen, asize)
        qdata = np.array(qdata).T

        # size distribution
        phi = arange ** self.pexp

        # normalizing mass
        cnorm = 4. / 3. * np.pi * self.rho_bulk * (self.amax ** pexp4
                                                   - self.amin ** pexp4) / pexp4
        # store inverse
        icnorm = 1e0 / cnorm

        kappa = []
        # loop on asize to compute integral
        for q in qdata:
            kappa.append(np.pi * np.trapz(q * arange ** 2 * phi, arange) * icnorm)

        # copy to attribute
        self.data["wlen"] = self.optical.data["wlen"]
        self.data["freq"] = self.optical.data["freq"]
        self.data["kappa"] = kappa

    # ***********************
    # compute opacity (cross section per dust mass), cm2/g
    # Assumes MNR dust size distribution with power-law exponent pexp,
    # in the range amin, amax, rho_bulk bulk density. Integral on dust
    # distribution uses ngrid points. The coating layer has thickness
    # aleyer
    def compute_kappa_coating(self, verbose=1):
        import numpy as np
        import sys

        # power-law exponent + 4 for integrals
        pexp4 = self.pexp + 4.

        # size range of grain core, cm
        arange = np.logspace(np.log10(self.amin), np.log10(self.amax), self.ngrid)

        if verbose > 0:
            print "Computing opacity with coating..."

        if self.alayer is None:
            sys.exit("ERROR: for composite materials you should set dust.alayer")

        qdata = []
        # loop on grain sizes
        for ii, asize in enumerate(arange):
            if verbose > 0:
                print round(ii * 1e2 / (self.ngrid - 1), 1), "%"
            asize_coat = asize + self.alayer
            self.optical.compute_q([asize, asize_coat])
            qabs = self.optical.data["qabs"]
            qdata.append(qabs)
        # transpose from (asize, wlen) to (wlen, asize)
        qdata = np.array(qdata).T

        # size range of grain core+mantle, cm
        arange_full = np.logspace(np.log10(self.amin + self.alayer),
                                  np.log10(self.amax + self.alayer), self.ngrid)

        # size distribution including mantle
        phi = arange_full ** self.pexp

        # normalizing mass
        cnorm = 4. / 3. * np.pi * self.rho_bulk * ((self.amax + self.alayer) ** pexp4
                                                   - (self.amin + self.alayer) ** pexp4) / pexp4
        icnorm = 1e0 / cnorm

        kappa = []
        # loop on wlen to compute integral on the size distribution
        for q in qdata:
            kappa.append(np.pi * np.trapz(q * arange_full ** 2 * phi, arange_full) * icnorm)

        # copy to attribute
        self.data["wlen"] = self.optical.data["wlen"]
        self.data["freq"] = self.optical.data["freq"]
        self.data["kappa"] = kappa

    # ****************
    def reverse_kappa(self):
        import sys
        import numpy as np
        from material import Material
        from optical import Optical
        nstep = 10
        ref_max = np.log10(2e0)
        ref_min = np.log10(1e-2)
        fout = open("surf.dat", "w")
        for ii, wlen in enumerate(self.data["wlen"][::100]):
            if wlen < 0.1:
                continue
            print ii
            kappa = self.data["kappa"][ii]
            for real_m in np.logspace(ref_min, ref_max, nstep):
                xmin = 1e1**ref_min
                xmax = 1e1**ref_max

                def get_kappa(rm, im_m):
                    mat = Material()
                    mat.add_value(wlen, complex(rm, im_m))
                    opt = Optical(mat)
                    opt.compute_kappa(verbose=0)
                    return opt.dust.data["kappa"][0]

                fmin = kappa - get_kappa(real_m, xmin)
                fmax = kappa - get_kappa(real_m, xmax)
                if fmin*fmax > 0e0:
                    print fmin, fmax, xmin, xmax
                    print("ERROR: same signs in bisection")
                    continue

                xmed = None
                while abs(xmax - xmin) > 1e-6:
                    xmed = (xmin + xmax) / 2e0
                    fmed = kappa - get_kappa(real_m, xmed)
                    if fmed*fmax < 0e0:
                        xmin = xmed
                        fmin = fmed
                    else:
                        xmax = xmed
                        fmax = fmed

                fout.write(str(wlen) + " " + str(real_m) + " " + str(xmed) + "\n")
                fout.flush()
            fout.write("\n")
