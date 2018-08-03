import math
import numpy as np
import re
import argparse


# TODO: empty prefix to the hasattr method
class PlotIgor(object):
    def __init__(self, infile, outfile=None, prefix=None):
        self.infile = infile
        self.outfile = outfile or (str(infile) + ".itx")
        self.prefix = prefix or ""
        self.wave = None
        return

    def read_bands(self):
        with open(self.infile, "r") as bandfile:
            kpts = []
            band = []
            kpath = []
            highsym = []

            index = bandfile.readline().split()
            numbands = int(index[2].strip(","))
            numkpts = int(index[4].strip())
            bandperline = int(np.ceil(numbands / 10))

            lines = bandfile.readlines()
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

            kpath = np.reshape(np.array(kpath, dtype='d'), (len(kpath), 1))
            highsym = np.reshape(np.array(highsym, dtype='d'), (len(highsym), 1))

            dic = {"kpts": kpts,
                   "kpath": kpath,
                   "highsym": highsym,
                   "band": band}

            self.wave = dic
            return

    def write_bands(self, plot=True, fermi=0.0, shift=False, guide=False):
        layout_preset = ("X DefaultFont/U \"Times New Roman\"\n"
                         "X ModifyGraph width=255.118,height=340.157\n"
                         "X ModifyGraph marker=19\n"
                         "X ModifyGraph lSize=1.5\n"
                         "X ModifyGraph tick(left)=2,tick(bottom)=3,noLabel(bottom)=2\n"
                         "X ModifyGraph mirror=1\n"
                         "X ModifyGraph zero(left)=8\n"
                         "X ModifyGraph fSize=28\n"
                         "X ModifyGraph lblMargin(left)=15,lblMargin(bottom)=10\n"
                         "X ModifyGraph standoff=0\n"
                         "X ModifyGraph axThick=1.5\n"
                         "X ModifyGraph axisOnTop=1\n"
                         "X Label left \"\Z28 Energy (eV)\"\n"
                         "X ModifyGraph zero(bottom)=0;DelayUpdate\n"
                         "X SetAxis left -3,3\n"
                         "X ModifyGraph zeroThick(left)=2.5\n"
                         )

        if self.prefix != "":
            waveprefix = str(self.prefix) + "_"
        else:
            waveprefix = input("Please type the system name : ") + "_"

        guide_preset = ("X AppendToGraph " + waveprefix + "guide_y1 " + waveprefix + "guide_y2 vs " +
                        waveprefix + "k_highsym\n"
                        "X ModifyGraph mode(" + waveprefix + "guide_y1)=1,rgb(" + waveprefix + "guide_y1)=(0,0,0)\n"
                        "X ModifyGraph mode(" + waveprefix + "guide_y2)=1,rgb(" + waveprefix + "guide_y2)=(0,0,0)\n"
                        "X SetAxis left -3,3"
                        )

        if shift is True:
            vbm = -10.0
            for x in self.wave["band"]:
                for y in x:
                    if (y <= fermi) and (vbm <= y):
                        vbm = y
            self.wave["band"] -= vbm
        else:
            self.wave["band"] -= fermi

        with open(self.outfile, "w") as out:
            # tmp = []
            out.write("IGOR\n")
            out.write("WAVES/D")
            out.write(" %s%s" % (waveprefix, "kpath"))
            for i in range(np.shape(self.wave["band"])[1]):
                out.write(" %s%s_%s" % (waveprefix, "band", i))
                # tmp.append("%s%s_%s" % (waveprefix, "band", i))
            out.write("\n")
            out.write("BEGIN\n")
            for i in range(len(self.wave["band"])):
                out.write(" %s" % self.wave["kpath"][i][0])
                for values in self.wave["band"][i]:
                    out.write(" %s" % values)
                out.write("\n")
            out.write("END\n")

            out.write("WAVES/D")
            out.write(" %s%s %s%s %s%s\n" % (waveprefix, "k_highsym", waveprefix, "guide_y1", waveprefix, "guide_y2"))
            out.write("BEGIN\n")
            for values in self.wave["highsym"]:
                out.write(" %s -30.00 30.00\n" % values[0])
            out.write("END\n")

            if plot is True:
                out.write("X Display %s%s_0 vs %s%s as \"%s%s\"\n" % (waveprefix, "band", waveprefix, "kpath",
                                                                      waveprefix, "band"))
                for i in range(np.shape(self.wave["band"])[1]):
                    if i == 0:
                        pass
                    else:
                        out.write("X AppendToGraph %s%s_%s vs %s%s\n" % (waveprefix, "band", i, waveprefix, "kpath"))
                out.write(layout_preset)

            if guide is True:
                out.write(guide_preset)

        return

    def read_dos(self):
        with open(self.infile, "r") as dosfile:
            egrid = []
            dos = []

            index = dosfile.readline().split()
            efermi = float(index[-2].strip())

            lines = dosfile.readlines()
            for x in lines:
                egrid.append(x.split()[0])
                dos.append(x.split()[1:])

        dic = {"egrid": np.array(egrid, dtype='d'),
               "dos": np.array(dos, dtype='d'),
               "efermi": efermi
               }

        self.wave = dic
        return

    def read_pdos(self):
        return

    def write_dos(self, plot=True, fermi=0.0):
        layout_preset = ("X DefaultFont/U \"Times New Roman\"\n"
                         "X ModifyGraph width=340.157,height=255.118\n"
                         "X ModifyGraph marker=19\n"
                         "X ModifyGraph lSize=1.5\n"
                         "X ModifyGraph tick(left)=2,tick(bottom)=3,noLabel(bottom)=2\n"
                         "X ModifyGraph mirror=1\n"
                         "X ModifyGraph zero(left)=8\n"
                         "X ModifyGraph fSize=28\n"
                         "X ModifyGraph lblMargin(left)=15,lblMargin(bottom)=10\n"
                         "X ModifyGraph standoff=0\n"
                         "X ModifyGraph axThick=1.5\n"
                         "X ModifyGraph axisOnTop=1\n"
                         "X Label bottom \"\Z28 Energy (eV)\"\n"
                         "X Label left \"\Z28 DOS (arb. unit)\"\n"
                         "X ModifyGraph zero(bottom)=0;DelayUpdate\n"
                         "X SetAxis bottom -3,3\n"
                         "X ModifyGraph zeroThick(left)=2.5\n"
                         )

        if self.prefix != "":
            waveprefix = str(self.prefix) + "_"
        else:
            waveprefix = input("Please type the system name : ") + "_"

        if np.shape(self.wave["dos"])[1] == 2:
            spin = False
        else:
            spin = True

        self.wave["Egrid"] -= fermi

        with open(self.outfile, "w") as out:
            out.write("IGOR\n")
            out.write("WAVES/D")
            out.write(" %s%s" % (waveprefix, "Egrid"))
            if spin is False:
                out.write(" %s%s" % (waveprefix, "tdos"))
                out.write(" %s%s" % (waveprefix, "intdos"))
            elif spin is True:
                out.write(" %s%s" % (waveprefix, "tdos_up"))
                out.write(" %s%s" % (waveprefix, "tdos_dw"))
                out.write(" %s%s" % (waveprefix, "intdos"))
            out.write("\n")
            out.write("BEGIN\n")
            for i in range(len(self.wave["dos"])):
                out.write(" %s" % self.wave["Egrid"][i][0])
                for values in self.wave["dos"][i]:
                    out.write(" %s" % values)
                out.write("\n")
            out.write("END\n")

            if plot is True:
                out.write("X Display %s%s vs %s%s as \"%s%s\"\n" %
                          (waveprefix, "tdos", waveprefix, "Egrid", waveprefix, "tdos"))
                out.write(layout_preset)

        return
