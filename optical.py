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
    def add_plot_kappa(self, fname, postfix="", linestyle="-", xlim=None, ylim=None, color=None,
                       linewidth=1, legend_columns=1):
        self.dust.add_plot_kappa(fname, postfix=postfix, linestyle=linestyle,
                                 xlim=xlim, ylim=ylim, color=color, linewidth=linewidth,
                                 legend_columns=legend_columns)

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
        for row in open(fname):
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
    def add_plot_q(self, fname, ptype="loglog", linestyle="-", postfix="", what=None,
                   colors=None):
        import matplotlib.pyplot as plt

        if colors is None:
            colors = ["tab:blue", "tab:orange", "tab:green"]

        if what is None:
            what = ["abs", "sca"]

        print("Plotting Q* to " + fname + "...")
        if type(ptype) is str:
            ptype = eval("plt." + ptype)

        for ii, w in enumerate(what):
            ptype(self.data["wlen"], self.data["q" + w], label="$Q_{" + w + "}$" + postfix,
                  linestyle=linestyle, color=colors[ii])

        plt.xlabel("$\\lambda$ / $\\mu$m")
        plt.ylabel("$Q_x$")
        plt.legend(loc="best")
        plt.savefig(fname)

    # ********************
    # plot refraction index
    def plot_ref_index(self, fname, ptype="loglog", linestyle="-", marker=None):
        import matplotlib.pyplot as plt
        plt.clf()
        self.add_plot_ref_index(fname, ptype=ptype, linestyle=linestyle, marker=marker)

    # ********************
    # plot refraction index
    def add_plot_ref_index(self, fname, ptype="loglog", linestyle="-", marker=None):
        import matplotlib.pyplot as plt

        print("Plotting refractive index to " + fname)

        if type(ptype) is str:
            ptype = eval("plt." + ptype)

        for material in self.materials:
            ptype(material.data["wlen"], material.data["real_m"], label="Re($m$)",
                  linestyle=linestyle, marker=marker)
            ptype(material.data["wlen"], material.data["im_m"], label="Im($m$)",
                  linestyle=linestyle, marker=marker)

        plt.xlabel("$\\lambda$ / $\\mu$m")
        plt.ylabel("$m$")
        plt.legend(loc="best")
        plt.savefig(fname)

    # ******************
    # save refrative index to file(s)
    # note: multiple material saved to multiple files
    def save_refractive_index(self, fname):
        import sys

        # convert finame to list if not a list
        if type(fname) is not list:
            fname = [fname]

        # check number of filenames
        if len(self.materials) != len(fname):
            sys.exit("ERROR: number of filenames must be same as materials!")

        # loop on materials to save refractive index to file
        for ii, material in enumerate(self.materials):
            material.save_refractive_index(fname[ii])

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
        data_q = {lab: np.array(x) for lab, x in data_q.items()}

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
            sys.exit("ERROR: the radius of the mantle cannot be smaller or equal "
                     "to the core radius!")

        # local copy of data to work with
        data = materials[0].data
        data_coat = materials[1].data

        # get interpolating function with optical properties
        # of the mantle material
        f_coat_real_m = interp1d(data_coat["wlen"], data_coat["real_m"])
        f_coat_im_m = interp1d(data_coat["wlen"], data_coat["im_m"])

        f_core_real_m = interp1d(data["wlen"], data["real_m"])
        f_core_im_m = interp1d(data["wlen"], data["im_m"])

        # get optical properties range of the optical properties
        # of the coating material
        wlen_min = max(data_coat["wlen"][0], data["wlen"][0])
        wlen_max = min(data_coat["wlen"][-1], data["wlen"][-1])

        # initialize data dictionary
        data_q_coat = {"wlen": [], "qabs": [], "qsca": [], "qbak": []}
        wlen_all = sorted(np.concatenate((data["wlen"], data_coat["wlen"])))
        # loop on wavelengths to compute Qx
        for ii, wlen in enumerate(wlen_all):
            # cannot calculate outside ranges
            if wlen < wlen_min or wlen > wlen_max:
                continue
            # prepare complex refractive index
            # ref_core = complex(data["real_m"][ii], data["im_m"][ii])
            ref_core = complex(f_core_real_m(wlen), f_core_im_m(wlen))
            ref_coat = complex(f_coat_real_m(wlen), f_coat_im_m(wlen))
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
        data_q_coat = {lab: np.array(x) for lab, x in data_q_coat.items()}
        self.data = data_q_coat

    # ****************************
    # add spherical impurities with given volume filling factor.
    # optical prpoerties of the impurities are taken from the material in the optical_impurity
    # object
    def add_impurity(self, optical_impurity_list, volume_filling_fraction_list):
        import sys

        for optical_impurity in optical_impurity_list:
            # rise error if there are more materials in the optical object
            if len(optical_impurity.materials) != 1:
                sys.exit("ERROR: impurity optical must have only one material!")

        # add optical impurities to the materials
        for material in self.materials:
            material.add_impurity([x.materials[0] for x in optical_impurity_list],
                                  volume_filling_fraction_list)

    # ****************************
    def extrapolate(self, wmax, nlast=45):

        for ii, material in enumerate(self.materials):
            self.materials[ii].extrapolate(wmax, nlast=nlast)
