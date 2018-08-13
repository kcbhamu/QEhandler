from qwutils.generalutils import Unitconverter
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

    @staticmethod
    def layout_preset(plottype, option=None):
        if plottype == "band":
            preset = ("X DefaultFont/U \"Times New Roman\"\n"
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
                      "X ModifyGraph zeroThick(left)=1.5\n"
                      )
        elif plottype == "dos":
            preset = ("X DefaultFont/U \"Times New Roman\"\n"
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
                      )

        elif plottype == "wf":
            preset = ("X DefaultFont/U \"Times New Roman\"\n"
                      "X ModifyGraph width=340.157,height=255.118\n"
                      "X ModifyGraph marker=19\n"
                      "X ModifyGraph lSize=1.5\n"
                      "X ModifyGraph tick=2\n"
                      "X ModifyGraph mirror=1\n"
                      "X ModifyGraph fSize=28\n"
                      "X ModifyGraph lblMargin(left)=15,lblMargin(bottom)=10\n"
                      "X ModifyGraph standoff=0\n"
                      "X ModifyGraph axThick=1.5\n"
                      "X ModifyGraph axisOnTop=1\n"
                      "X Label bottom \"\Z28 Distance (\{num2char(129)})\"\n"
                      "X Label left \"\Z28 Energy (eV)\"\n"
                      "X ModifyGraph zero(bottom)=0;DelayUpdate\n"
                      )

        elif plottype == "diel":
            preset = ("X DefaultFont/U \"Times New Roman\"\n"
                      "X ModifyGraph width=340.157,height=255.118\n"
                      "X ModifyGraph marker=19\n"
                      "X ModifyGraph lSize=1.5\n"
                      "X ModifyGraph tick=2\n"
                      "X ModifyGraph mirror=1\n"
                      "X ModifyGraph fSize=28\n"
                      "X ModifyGraph lblMargin(left)=15,lblMargin(bottom)=10\n"
                      "X ModifyGraph standoff=0\n"
                      "X ModifyGraph axThick=1.5\n"
                      "X ModifyGraph axisOnTop=1\n"
                      "X Label bottom \"\Z28Photon Energy (eV))\"\n"
                      "X Label left \"\Z28\F'Symbol'e\B%s\M\F'Times New Roman' (a.u.)\"\n"
                      "X ModifyGraph zero(bottom)=0;DelayUpdate\n"
                      "X ModifyGraph zero(left)=8,zeroThick(left)=1.5\n"
                      % option)

        return preset

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
                out.write(" %s%s_%s" % (waveprefix, "band", i + 1))
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
                out.write(self.layout_preset("band"))

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

        egrid = np.array(egrid, dtype='d')
        dos = np.array(dos, dtype='d')

        dic = {"dos": {"egrid": np.reshape(egrid, (1, len(egrid), 1)),
                       "dos": np.reshape(dos, (1, np.shape(dos)[0], np.shape(dos)[1])),
                       "efermi": efermi
                       }
               }
        self.wave = dic
        return

    def read_pdos(self, kproj=False):
        """
* The format for the collinear, spin-unpolarized case and the
  non-collinear, spin-orbit case is:
      E DOS(E) PDOS(E)
      ...

* The format for the collinear, spin-polarized case is:
      E DOSup(E) DOSdw(E)  PDOSup(E) PDOSdw(E)
      ...

* The format for the non-collinear, non spin-orbit case is:
      E DOS(E) PDOSup(E) PDOSdw(E)
      ...

In the collinear case and the non-collinear, non spin-orbit case
projected DOS are written to file "filpdos".pdos_atm#N(X)_wfc#M(l),
where N = atom number , X = atom symbol, M = wfc number, l=s,p,d,f
(one file per atomic wavefunction found in the pseudopotential file)

* The format for the collinear, spin-unpolarized case is:
      E LDOS(E) PDOS_1(E) ... PDOS_2l+1(E)
      ...
  where LDOS = \sum m=1,2l+1 PDOS_m(E)
  and PDOS_m(E) = projected DOS on atomic wfc with component m

* The format for the collinear, spin-polarized case and the
  non-collinear, non spin-orbit case is as above with
  two components for both  LDOS(E) and PDOS_m(E)

In the non-collinear, spin-orbit case (i.e. if there is at least one
fully relativistic pseudopotential) wavefunctions are projected
onto eigenstates of the total angular-momentum.
Projected DOS are written to file "filpdos".pdos_atm#N(X)_wfc#M(l_j),
where N = atom number , X = atom symbol, M = wfc number, l=s,p,d,f
and j is the value of the total angular momentum.
In this case the format is:
    E LDOS(E) PDOS_1(E) ... PDOS_2j+1(E)
    ...

If kresolveddos=.true., the k-point index is prepended
to the formats above, e.g. (collinear, spin-unpolarized case)
    ik E DOS(E) PDOS(E)

All DOS(E) are in states/eV plotted vs E in eV

Orbital Order

 Order of m-components for each l in the output:

    1, cos(phi), sin(phi), cos(2*phi), sin(2*phi), .., cos(l*phi), sin(l*phi)

where phi is the polar angle:x=r cos(theta)cos(phi), y=r cos(theta)sin(phi)
This is determined in file Modules/ylmr2.f90 that calculates spherical harmonics.

for l=1:
  1 pz     (m=0)
  2 px     (real combination of m=+/-1 with cosine)
  3 py     (real combination of m=+/-1 with sine)

for l=2:
  1 dz2    (m=0)
  2 dzx    (real combination of m=+/-1 with cosine)
  3 dzy    (real combination of m=+/-1 with sine)
  4 dx2-y2 (real combination of m=+/-2 with cosine)
  5 dxy    (real combination of m=+/-2 with sine)

        """

        with open(self.infile, "r") as pdosfile:
            egrid = []
            dos = []

            lines = pdosfile.readlines()
            for x in lines:
                if len(x) == 0:
                    pass
                elif kproj is False:
                    egrid.append(x.split()[0])
                    dos.append(x.split()[1:])
                elif kproj is True:
                    kp.append(x.split()[0])
                    egrid.append(x.split()[1])
                    dos.append(x.split()[2:])

        egrid = np.array(egrid, dtype='d')
        dos = np.array(dos, dtype='d')

        return

    def write_dos(self, plot=True, fermi=0.0):
        if self.prefix != "":
            waveprefix = str(self.prefix) + "_"
        else:
            waveprefix = input("Please type the system name : ") + "_"

        if np.shape(self.wave["dos"])[1] == 2:
            spin = False
        else:
            spin = True

        self.wave["egrid"] -= fermi

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
                out.write(" %s" % self.wave["egrid"][i][0])
                for values in self.wave["dos"][i]:
                    out.write(" %s" % values)
                out.write("\n")
            out.write("END\n")

            if plot is True:
                if spin is True:
                    out.write("X Display %s%s vs %s%s as \"%s%s\"\n" %
                              (waveprefix, "tdos_up", waveprefix, "Egrid", waveprefix, "tdos"))
                elif spin is False:
                    out.write("X Display %s%s vs %s%s as \"%s%s\"\n" %
                              (waveprefix, "tdos", waveprefix, "Egrid", waveprefix, "tdos"))
                out.write(self.layout_preset("dos"))

        return

    def read_wf(self):
        with open(self.infile, "r") as wffile:
            distance = []
            macro = []
            planar = []
            lines = wffile.readlines()

            for x in lines:
                distance.append(Unitconverter.unit_convert(float(x.split()[0]), "length", "bohr", "ang"))
                planar.append(Unitconverter.unit_convert(float(x.split()[1]), "energy", "Ry", "eV"))
                macro.append(Unitconverter.unit_convert(float(x.split()[2]), "energy", "Ry", "eV"))

            distance = np.array(distance, dtype='d')
            macro = np.array(macro, dtype='d')
            planar = np.array(planar, dtype='d')

            dic = {"distance": np.reshape(distance, (len(distance), 1)),
                   "macro": np.reshape(macro, (len(macro), 1)),
                   "planar": np.reshape(planar, (len(planar), 1)),
                   }
        self.wave = dic
        return

    def write_wf(self, plot=True):
        with open(self.outfile, "w") as out:

            if self.prefix != "":
                waveprefix = str(self.prefix) + "_"
            else:
                waveprefix = input("Please type the system name : ") + "_"

            out.write("IGOR\n")
            out.write("WAVES/D")
            out.write(" %s%s  %s%s  %s%s\n" % (waveprefix, "Distance", waveprefix, "Macroavg", waveprefix, "Planaravg"))
            out.write("BEGIN\n")
            for i in range(len(self.wave["distance"])):
                out.write(" %s" % self.wave["distance"][i][0])
                out.write(" %s" % self.wave["macro"][i][0])
                out.write(" %s" % self.wave["planar"][i][0])
                out.write("\n")
            out.write("END\n")

            if plot is True:
                out.write("X Display %s%s vs %s%s as \"%s%s\"\n" %
                          (waveprefix, "Macroavg", waveprefix, "Distance", waveprefix, "potential"))
                out.write("X AppendToGraph %s%s vs %s%s\n" % (waveprefix, "Planaravg", waveprefix, "Distance"))
                out.write(self.layout_preset("wf"))

        return

    def read_diel(self, real=None, imag=None, ieps=None, eels=None, direction=False):
        dic = {}

        def reader(infile):
            with open(infile, "r") as file:
                e = []
                x = []
                y = []
                z = []
                lines = file.readlines()

                for i in range(len(lines) - 2):
                    e.append(lines[i + 2].split()[0])
                    x.append(lines[i + 2].split()[1])
                    y.append(lines[i + 2].split()[2])
                    z.append(lines[i + 2].split()[3])

                e = np.array(e, dtype='d')
                x = np.array(x, dtype='d')
                y = np.array(y, dtype='d')
                z = np.array(z, dtype='d')
                avg = (x + y + z) / 3

            return e, avg, x, y, z

        if real is not None:
            tmp = reader(real)
            dic["e1_E"] = tmp[0]
            dic["e1"] = tmp[1]
            if direction is True:
                dic["e1_x"] = tmp[2]
                dic["e1_y"] = tmp[3]
                dic["e1_z"] = tmp[4]

        if imag is not None:
            tmp = reader(imag)
            dic["e2_E"] = tmp[0]
            dic["e2"] = tmp[1]
            if direction is True:
                dic["e2_x"] = tmp[2]
                dic["e2_y"] = tmp[3]
                dic["e2_z"] = tmp[4]

        if ieps is not None:
            tmp = reader(ieps)
            dic["diag_E"] = tmp[0]
            dic["diag"] = tmp[1]
            if direction is True:
                dic["diag_x"] = tmp[2]
                dic["diag_y"] = tmp[3]
                dic["diag_z"] = tmp[4]

        if eels is not None:
            tmp = reader(eels)
            dic["eels_E"] = tmp[0]
            dic["eels"] = tmp[1]
            if direction is True:
                dic["eels_x"] = tmp[2]
                dic["eels_y"] = tmp[3]
                dic["eels_z"] = tmp[4]

        self.wave = dic
        return

    def write_diel(self, plot=True):
        with open(self.outfile, "w") as out:
            wavename = []
            if self.prefix != "":
                waveprefix = str(self.prefix) + "_"
            else:
                waveprefix = input("Please type the system name : ") + "_"

            out.write("IGOR\n")
            out.write("WAVES/D")
            for x in self.wave.keys():
                out.write(" %s%s" % (waveprefix, x))
                wavename.append(x)
            out.write("\n")
            out.write("BEGIN\n")

            for i in range(len(self.wave[wavename[0]])):
                for x in wavename:
                    out.write(" %s" % self.wave[x][i])
                out.write("\n")
            out.write("END\n")

            if plot is True:
                out.write("X Display %s%s vs %s%s as \"%s%s\"\n" %
                          (waveprefix, "e1", waveprefix, "e1_E", waveprefix, "e1"))
                out.write(self.layout_preset("diel", 1))
                out.write("X Display %s%s vs %s%s as \"%s%s\"\n" %
                          (waveprefix, "e2", waveprefix, "e2_E", waveprefix, "e2"))
                out.write(self.layout_preset("diel", 2))
        return
