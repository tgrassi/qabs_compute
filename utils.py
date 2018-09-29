from material import Material
from optical import Optical


class QabsManager:

    # ******************
    def __init__(self):
        self.materials = dict()
        self.opticals = dict()

    # ********************
    def load_material(self, fname, labs=None, name=None, force=False):
        import sys

        # default labs input
        default_labs = ["wlen", "real_eps1", "im_eps", "real_m1", "im_m"]

        # use default is labs is not set
        if labs is None:
            labs = default_labs

        # generate name automatically
        if name is None:
            name = "material_" + str(len(self.materials))

        if name in self.materials and not force:
            print "ERROR: material " + name + " already present!"
            print "Change name or use force=True"
            sys.exit()

        mat = Material(fname, labs, name)
        opt = Optical(mat)

        # add materials to dictionary
        self.materials[name] = mat
        self.opticals[name] = opt

        return opt

    # *****************
    # create new optical material from others
    def make_optical(self, opticals, name=None):
        import sys

        # check input
        if len(opticals) != 2:
            sys.exit("ERROR: you need 2 materials to create a composite one!")

        # generate a random name
        if name is None:
            name = "optical_" + str(len(self.opticals))

        # loop on opticals to concatenate materials
        materials = []
        for opt in opticals:
            materials += opt.materials

        # store optical to attribute dictionary
        opt = Optical(materials)
        self.opticals[name] = opt
        return opt

    # ******************
    # write a database report
    def report(self):
        # loop on standard data and coating
        for name, material in self.materials.iteritems():
            data = material.data
            # if missing data skip
            if "wlen" not in data:
                continue
            # small table with data
            print "*********"
            print "Name:", name
            print "File:", data["fname"]
            print "Number of points:", len(data["wlen"])
            print "Range wavelength (micron): %e, %e" % (min(data["wlen"]), max(data["wlen"]))
            print "Range frequency (Hz): %e, %e" % (min(data["freq"]), max(data["freq"]))

    # ************************
    # perform a benchmark
    def benchmark(self, fname="benchmark_q.png"):

        # load material
        opt = self.load_material("data/eps_Sil.dat", labs=["wlen", "real_eps1", "im_eps",
                                                           "real_m1", "im_m"])

        # compute Qabs for given size, cm
        opt.compute_q(1e-7)
        # plot computed
        opt.plot_q(fname)

        # load Qabs from file to compare
        opt.load_q("data/Sil_21_1e3.dat")
        # over-plot loaded
        opt.add_plot_q(fname, linestyle="--")

    # *******************
    @staticmethod
    def clear_plots():
        import matplotlib.pyplot as plt
        plt.clf()