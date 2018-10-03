from dust import Dust


class Optical:

    # *****************
    # class constructor with materials
    def __init__(self, materials):

        if type(materials) is not list:
            materials = [materials]

        self.materials = materials
        self.data = dict()

        self.dust = Dust(self)

    # **************************
    # plot opacity, postfix will be added in the legend
    def plot_kappa(self, fname, postfix="", linestyle="-"):
        self.dust.plot_kappa(fname, postfix=postfix, linestyle=linestyle)

    # **************************
    # add opacity plot to current plt session, postfix will be added in the legend
    def add_plot_kappa(self, fname, postfix="", linestyle="-"):
        self.dust.add_plot_kappa(fname, postfix=postfix, linestyle=linestyle)

    # ******************
    # save opacity to file
    def save_kappa(self, fname):
        self.dust.save_kappa(fname)

    # ******************
    # load opacity from file
    def load_kappa(self, fname):
        self.dust.load_kappa(fname)

    # **************************
    # compute opacity
    def compute_kappa(self, verbose=1):
        self.dust.compute_kappa(verbose)

    # *******************
    # load optical properties from file
    def load_q(self, fname, labs=None):

        # default labels
        if labs is None:
            labs = ["wlen", "qabs", "qsca", "gcos"]

        self.data = {ll: [] for ll in labs}
        for row in open(fname, "rb"):
            srow = row.strip()
            if srow.startswith("#"):
                continue
            if srow == "":
                continue
            arow = [float(x) for x in srow.split(" ") if x != ""]
            for ii, ll in enumerate(labs):
                self.data[ll].append(arow[ii])

    # **************************
    # plot optical properties to file, postfix is for the legend
    def plot_q(self, fname, ptype="loglog", linestyle="-", postfix=""):
        import matplotlib.pyplot as plt
        plt.clf()
        self.add_plot_q(fname, ptype=ptype, linestyle=linestyle, postfix=postfix)

    # **************************
    # add optical properties plot
    def add_plot_q(self, fname, ptype="loglog", linestyle="-", postfix=""):
        import matplotlib.pyplot as plt

        print "Plotting Q* to " + fname + "..."
        if type(ptype) is str:
            ptype = eval("plt." + ptype)

        ptype(self.data["wlen"], self.data["qabs"], label="$Q_{abs}$" + postfix,
              linestyle=linestyle)
        ptype(self.data["wlen"], self.data["qsca"], label="$Q_{sca}$" + postfix,
              linestyle=linestyle)

        plt.xlabel("$\\lambda$ / $\\mu$m")
        plt.ylabel("$Q_x$")
        plt.legend(loc="best")
        plt.savefig(fname)

    # ********************
    # plot refraction index
    def plot_ref_index(self, fname, ptype="loglog"):
        import matplotlib.pyplot as plt
        plt.clf()
        self.add_plot_ref_index(fname, ptype=ptype)

    # ********************
    # plot refraction index
    def add_plot_ref_index(self, fname, ptype="loglog"):
        import matplotlib.pyplot as plt

        print "Plotting refractive index to " + fname

        if type(ptype) is str:
            ptype = eval("plt." + ptype)

        for material in self.materials:
            ptype(material.data["wlen"], material.data["real_m"], label="Re($m$)")
            ptype(material.data["wlen"], material.data["im_m"], label="Im($m$)")

        plt.xlabel("$\\lambda$ / $\\mu$m")
        plt.ylabel("$m$")
        plt.legend(loc="best")
        plt.savefig(fname)

    # *********************
    # compute Q* for a single size
    def compute_q(self, asizes):
        import sys

        if type(asizes) is not list:
            asizes = [asizes]

        # local variable
        materials = self.materials

        # check number of sizes
        if len(materials) != len(asizes):
            sys.exit("ERROR: number of sizes does not match number of materials!")

        # check what subroutine to use
        if len(materials) == 1:
            self.compute_q_bare(materials[0], asizes[0])
        elif len(materials) == 2:
            self.compute_q_coating(materials, asizes)
        else:
            sys.exit("ERROR: you cannot combine more than two materials!")

    # ************************
    # compute qabs and qsca
    # data is the data structure loaded by load_eps
    def compute_q_bare(self, material, asize, nang=100):
        from bhmie import bhmie
        from utility import micron_to_hz
        import numpy as np

        # local copy of data to work with
        data = material.data

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
        data_q["freq"] = micron_to_hz(data_q["wlen"])

        # convert to numy arrays
        data_q = {lab: np.array(x) for lab, x in data_q.iteritems()}

        # copy data to attribute
        self.data = data_q

    # ***********
    # compute qabs, qsca, qback using coating
    # asizes is [core radius, mantle radius]
    # radius of the mantle
    def compute_q_coating(self, materials, asizes):
        from bhcoat import bhcoat_ph
        from scipy.interpolate import interp1d
        import numpy as np
        from utility import micron_to_hz
        import sys

        # check input radii
        if asizes[1] <= asizes[0]:
            sys.exit("ERROR: the radius of the mantle cannot be smaller or equal the core!")

        # local copy of data to work with
        data = materials[0].data
        data_coat = materials[1].data

        # get interpolating function with optical properties
        # of the mantle material
        f_real_m = interp1d(data_coat["wlen"], data_coat["real_m"])
        f_im_m = interp1d(data_coat["wlen"], data_coat["im_m"])

        # get optical properties range of the optical properties
        # of the coating material
        wlen_min = data_coat["wlen"][0]
        wlen_max = data_coat["wlen"][-1]

        # initialize data dictionary
        data_q_coat = {"wlen": [], "qabs": [], "qsca": [], "qbak": []}
        # loop on wavelengths to compute Qx
        for ii, wlen in enumerate(data["wlen"]):
            # cannot calculate outside ranges
            if wlen < wlen_min or wlen > wlen_max:
                continue
            # prepare complex refractive index
            ref_core = complex(data["real_m"][ii], data["im_m"][ii])
            ref_coat = complex(f_real_m(wlen), f_im_m(wlen))
            # call function to compute Qx, note wlen converted to cm
            qext, qsca, qbak = bhcoat_ph(asizes[0], asizes[1], ref_core, ref_coat, wlen*1e-4)

            # store output and compute derived quantities
            data_q_coat["qabs"].append(qext-qsca)
            data_q_coat["qsca"].append(qsca)
            data_q_coat["qbak"].append(qbak)
            data_q_coat["wlen"].append(wlen)

        # convert to Hz
        data_q_coat["freq"] = micron_to_hz(data_q_coat["wlen"])

        # convert to numy arrays
        data_q_coat = {lab: np.array(x) for lab, x in data_q_coat.iteritems()}
        self.data = data_q_coat
