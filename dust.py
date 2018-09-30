

class Dust:

    # ********************
    # dust class computes kappa and plot
    def __init__(self, optical=None, amin=5e-7, amax=2.5e-5, pexp=-3.5, rho_bulk=2.9,
                 ngrid=30, alayer=None):
        self.optical = optical  # dust material
        self.amin = amin  # min grain size, cm
        self.amax = amax  # max grain size, cm
        self.pexp = pexp  # MRN power-law exponent
        self.rho_bulk = rho_bulk  # bulk density, g/cm3
        self.ngrid = ngrid  # number of grid points for grain size (log-spaced)
        self.alayer = alayer  # mantle layer thickness (if any)
        self.data = dict()  # computed data

    # ******************
    # print report to screen
    def report(self):
        print "**************"
        print "Dust report"
        print "Material(s):", [x.name for x in self.optical.materials]
        print "min size, cm: %e" % self.amin
        print "max size, cm: %e" % self.amax
        print "MNR expoent: %f" % self.pexp
        print "bulk density, g/cm3: %e" % self.rho_bulk
        print "bins:", self.ngrid
        print "coating layer (if any), cm: %e" % self.alayer
        print "**************"

    # ******************
    # plot opacity to file
    def plot_kappa(self, fname, postfix=""):
        import matplotlib.pyplot as plt
        plt.clf()
        self.add_plot_kappa(fname, postfix=postfix)

    # *******************
    # add a plot to the current opacity plot
    def add_plot_kappa(self, fname, postfix=""):
        import matplotlib.pyplot as plt

        print "Plotting kappa to " + fname + "..."

        plt.loglog(self.data["wlen"], self.data["kappa"], label="$\\kappa$" + postfix)
        plt.xlabel("$\\lambda$ / $\\mu$m")
        plt.ylabel("$\\kappa$ / [cm$^2$ g$^{-1}$]")
        plt.legend(loc="best")
        plt.savefig(fname)

    # ******************
    # save opacity data to file
    def save_kappa(self, fname):

        fout = open(fname, "w")
        fout.write("# wlen(micron), opacity(cm2/g)\n")
        for ii, wlen in enumerate(self.data["wlen"]):
            fout.write(str(wlen) + " " + str(self.data["kappa"][ii]) + "\n")
        fout.close()

    # *********************
    # load data from file
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
    # compute opacity for the current dust
    def compute_kappa(self, verbose=1):
        import sys

        # check number of materials
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
        phi = arange_full**self.pexp

        # normalizing mass (without mantle, since mantle is supposed to be added to
        # the original dust distribution)
        cnorm = 4. / 3. * np.pi * self.rho_bulk * (self.amax**pexp4 - self.amin**pexp4) / pexp4
        icnorm = 1e0 / cnorm

        kappa = []
        # loop on wlen to compute integral on the size distribution
        for q in qdata:
            kappa.append(np.pi * np.trapz(q * arange_full**2 * phi, arange_full) * icnorm)

        # copy data to attribute variable
        self.data["wlen"] = self.optical.data["wlen"]
        self.data["freq"] = self.optical.data["freq"]
        self.data["kappa"] = kappa
