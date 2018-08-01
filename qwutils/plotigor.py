import math
import numpy as np
import re
import argparse


class PlotIgor(object):
    def __init__(self, infile, outfile=None):
        self.infile = infile
        self.outfile = outfile or (str(infile) + ".itx")
        self.wave = None
        return

    def readbands(self):
        with open(self.infile, "r") as file:
            kpts = []
            band = []
            kpath = []
            highsym = []

            index = file.readline().split()
            numbands = int(index[2].strip(","))
            numkpts = int(index[4].strip())
            bandperline = int(np.ceil(numbands / 10))

            lines = file.readlines()
            for i in range(numkpts):
                kpts.append(lines[i * (bandperline + 1)].split())

            for i in range(numkpts):
                tmp = []
                for j in range(bandperline):
                    tmp.append(lines[(i * (bandperline + 1)) + (j + 1)])
                band.append(" ".join(tmp).split())

            kpts = np.array(kpts, dtype='d')
            band = np.array(band, dtype='d')

            for i in range(len(kpts)):
                if i == 0:
                    kpath.append("0.0")
                    highsym.append("0.0")
                else:
                    kpath.append(float(kpath[-1]) + np.linalg.norm(kpts[i] - kpts[i - 1]))
                    if np.linalg.norm(kpts[i] - kpts[i - 1]) == 0:
                        highsym.append(float(kpath[-1]) + np.linalg.norm(kpts[i] - kpts[i - 1]))

            kpath = np.array(kpath, dtype='d')
            highsym = np.array(highsym, dtype='d')

            dic = {"kpts": kpts,
                   "kpath": kpath,
                   "highsym": highsym,
                   "band": band}

            self.wave = dic
            return

    def write_itx(self, plot=True):
        layout_preset = ("X DefaultFont/U \"Times New Roman\"\n"
                         "X ModifyGraph width=340.157,height=226.772\n"
                         "X ModifyGraph marker=19\n"
                         "X ModifyGraph lSize=1.5\n"
                         "X ModifyGraph tick=2\n"
                         "X ModifyGraph mirror=1\n"
                         "X ModifyGraph fSize=24\n"
                         "X ModifyGraph lblMargin(left)=15,lblMargin(bottom)=10\n"
                         "X ModifyGraph standoff=0\n"
                         "X ModifyGraph axThick=1.5\n"
                         "X ModifyGraph axisOnTop=1\n"
                         )

        with open(self.outfile, "w") as out:
            out.write("IGOR\n")
            out.write("WAVES/D")
            if self.header is not None:
                for x in self.header:
                    out.write(" " + str(x))
            else:
                for i in range(self.numitems):
                    out.write(" wave%s" % i)
            out.write("\n")
            out.write("BEGIN\n")

            for x in self.waves:
                for y in x:
                    out.write(" " + str(y))
                out.write("\n")
            out.write("END\n")

            if plot is True:
                if self.header is not None:
                    count = 0
                    for i in range(self.numitems):
                        if xindex == i:
                            pass
                        else:
                            if count == 0:
                                out.write("X Display %s vs %s\n" % (self.header[i], self.header[xindex]))
                                count += 1
                            else:
                                out.write("X AppendToGraph %s vs %s\n" % (self.header[i], self.header[xindex]))

                else:
                    count = 0
                    for i in range(self.numitems):
                        if xindex == i:
                            pass
                        else:
                            if count == 0:
                                out.write("X Display wave%s vs wave%s\n" % (i, xindex))
                                count += 1
                            else:
                                out.write("X AppendToGraph wave%s vs wave%s\n" % (i, xindex))

                out.write(layout_preset)

        return

    def plotband(self, plot2d=False):
        return
