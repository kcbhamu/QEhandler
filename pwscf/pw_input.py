from __future__ import unicode_literals, print_function
from collections import OrderedDict
from utils.generalutils import ScientificConstants, Unitconverter
import numpy as np
import yaml
import spglib


# TODO: writing all namelist tags at the end of class calling to prevent overlapping writing methods
class PWin(object):

    def __init__(self, file=None, cellparam=None, atompos=None, aspecies=None, kpoints=None, control=None, system=None,
                 electrons=None, ions=None, cell=None, constr=None, occup=None, aforces=None):
        self.cellparam = cellparam or {}
        self.atompos = atompos or {}
        self.aspecies = aspecies or {}
        self.kpoints = kpoints or {}
        self.control = control or {}
        self.system = system or {}
        self.electrons = electrons or {}
        self.ions = ions or {}
        self.cell = cell or {}
        self.constr = constr or {}
        self.occup = occup or {}
        self.aforces = aforces or {}
        self.infile = file

        with open("input_tags.yaml", 'r') as y:
            self.taglist = ordered_load(y, yaml.SafeLoader)[0]

    # TODO: handling namelists under "others" list
    def read_file(self):
        with open(self.infile, "r") as f:
            namelist = []
            others = ["ATOMIC_SPECIES", "CONSTRAINTS", "OCCUPATIONS", "ATOMIC_FORCES"]
            dic = OrderedDict()
            lines = f.readlines()
            for x in lines:
                line = x.strip("\n").strip(",")
                if "&" in line:
                    namelist.append(line)
                    dic[line] = OrderedDict()
                elif line in others:
                    namelist.append(line)
                    dic[line] = OrderedDict()
                elif "ATOMIC_POSITIONS" in line:
                    namelist.append(line.split()[0])
                    atom_pos = []
                    dic[namelist[-1]] = OrderedDict()
                    dic[namelist[-1]]["unit"] = line.split()[1]
                elif "CELL_PARAMETERS" in line:
                    namelist.append(line.split()[0])
                    cell_param = []
                    dic[namelist[-1]] = OrderedDict()
                    dic[namelist[-1]]["unit"] = line.split()[1]
                elif "K_POINTS" in line:
                    kp_list = []
                    namelist.append(line.split()[0])
                    dic[namelist[-1]] = OrderedDict()
                    dic[namelist[-1]]["unit"] = line.split()[1]
                elif line == "":
                    pass
                elif line == "/":
                    pass
                else:
                    if line[0] in "!#":
                        pass
                    elif namelist[-1] == "ATOMIC_SPECIES":
                        dic[namelist[-1]][pseudo_parser(line)[0]] = pseudo_parser(line)[1]
                    elif namelist[-1] == "CELL_PARAMETERS":
                        cell_param.append(line.split())
                    elif namelist[-1] == "K_POINTS":
                        kp_list.append(line.split())
                    elif namelist[-1] == "ATOMIC_POSITIONS":
                        atom_pos.append([len(atom_pos) + 1, coord_parser(line)])
                    else:
                        dic[namelist[-1]][tag_parser(line)[0]] = tag_parser(line)[1]

            if "&CONTROL" in namelist:
                self.control = dic["&CONTROL"]
            if "&SYSTEM" in namelist:
                self.system = dic["&SYSTEM"]
            if "&ELECTRONS" in namelist:
                self.electrons = dic["&ELECTRONS"]
            if "&IONS" in namelist:
                self.ions = dic["&IONS"]
            if "&CELL" in namelist:
                self.cell = dic["&CELL"]

            if "ATOMIC_SPECIES" in namelist:
                self.aspecies = dic["ATOMIC_SPECIES"]

            if "K_POINTS" in namelist:
                if dic["K_POINTS"]["unit"] == "automatic":
                    dic["K_POINTS"]["grid"] = [kp_list[0][0], kp_list[0][1], kp_list[0][2]]
                    dic["K_POINTS"]["shift"] = [kp_list[0][3], kp_list[0][4], kp_list[0][5]]
                elif dic["K_POINTS"]["unit"] == "gamma":
                    pass
                else:
                    dic["K_POINTS"]["numkp"] = kp_list.pop(0)
                    dic["K_POINTS"]["points"] = np.array(kp_list, dtype='d')
                self.kpoints = dic["K_POINTS"]

            if "CELL_PARAMETERS" in namelist:
                dic["CELL_PARAMETERS"]["vector"] = np.array(cell_param, dtype='d')
                self.cellparam = dic["CELL_PARAMETERS"]

            if "ATOMIC_POSITIONS" in namelist:
                elements = []
                coords = []
                dynamics = []
                for x in atom_pos:
                    elements.append(x[1][0])
                    coords.append(x[1][1])
                    dynamics.append(x[1][2])
                dic["ATOMIC_POSITIONS"]["elements"] = elements
                dic["ATOMIC_POSITIONS"]["coordinates"] = np.array(coords, dtype='d')
                dic["ATOMIC_POSITIONS"]["dynamics"] = np.array(dynamics, dtype='d')
                self.atompos = dic["ATOMIC_POSITIONS"]

        return

    def generate_pwin(self, tagvaluepair=None, pseudo=None, kpoints=None, cellparam=None, atomicpos=None):
        control = [("calculation", "scf"),
                   ("outdir", "./out/"),
                   ("etot_conv_thr", "1.0D-6"),
                   ("forc_conv_thr", "1.0D-5")
                   ]

        system = [("ibrav", "0"),
                  ("nat", ""),
                  ("ntyp", ""),
                  ("ecutwfc", "40.0"),
                  ]

        electrons = [("conv_thr", "1.0D-6")]

        self.control = OrderedDict(control)
        self.system = OrderedDict(system)
        self.electrons = OrderedDict(electrons)

        if tagvaluepair is not None:
            for i in range(len(tagvaluepair)):
                tag = value = None
                if i % 2 == 0:
                    tag = tagvaluepair[i]
                elif i % 2 == 0:
                    tag = tagvaluepair[i]
                self.write_tags(tag, value)

        # if pseudo is not None:

        return

    def write_tags(self, tag, value):
        for x in list(self.taglist.keys()):
            if tag in list(self.taglist[x].keys()):
                if x == "control":
                    self.control[tag] = value
                elif x == "system":
                    self.system[tag] = value
                elif x == "electrons":
                    self.electrons[tag] = value
                elif x == "ions":
                    self.ions[tag] = value
                elif x == "cell":
                    self.cell[tag] = value
            else:
                raise IOError("No such tag implemented in PWscf package!")

        return

    def remove_tags(self, tag, default=False):
        for x in list(self.taglist.keys()):
            if tag in list(self.taglist[x].keys()):
                if x == "control":
                    if default:
                        self.control[tag] = self.taglist[x][tag]['default']
                    else:
                        del self.control[tag]
                elif x == "system":
                    if default:
                        self.system[tag] = self.taglist[x][tag]['default']
                    else:
                        del self.system[tag]
                elif x == "electrons":
                    if default:
                        self.electrons[tag] = self.taglist[x][tag]['default']
                    else:
                        del self.electrons[tag]
                elif x == "ions":
                    if default:
                        self.ions[tag] = self.taglist[x][tag]['default']
                    else:
                        del self.ions[tag]
                elif x == "cell":
                    if default:
                        self.cell[tag] = self.taglist[x][tag]['default']
                    else:
                        del self.cell[tag]
            else:
                raise IOError("No such tag implemented in PWscf package!")
        return

    def add_pseudo(self, elem, mass=None, pseudo=None, pseudo_label=None):
        if mass is None:
            m = ScientificConstants().elementary_data(elem, "mass")
        else:
            m = mass

        if pseudo is None:
            if pseudo_label is None:
                pp = str(elem + "_dummy.upf")
            else:
                pp = str(elem + "_" + pseudo_label)
        else:
            pp = pseudo

        self.aspecies[elem] = {"mass": m, "pseudo": pp}

        return

    # TODO: test kpt file import
    def change_kpoints(self, unit="automatic", spacing=None, grid=None, shift=None, numkp=None, kpts=None, reduce=None):
        self.kpoints["unit"] = unit

        if unit == "gamma":
            pass
        elif unit == "automatic":
            self.kpoints["grid"] = grid
            self.kpoints["shift"] = shift
        elif unit == "spacing":
            self.kpoints["unit"] = "automatic"
            grida = int(np.ceil(np.linalg.norm(rec_vector(self.cellparam["vector"])[0]) * 2 * np.pi / spacing))
            gridb = int(np.ceil(np.linalg.norm(rec_vector(self.cellparam["vector"])[1]) * 2 * np.pi / spacing))
            gridc = int(np.ceil(np.linalg.norm(rec_vector(self.cellparam["vector"])[2]) * 2 * np.pi / spacing))

            if reduce is not None:
                if reduce == 0:
                    grida = 1
                elif reduce == 1:
                    gridb = 1
                elif reduce == 2:
                    gridc = 1

            self.kpoints["grid"] = "%s %s %s" % (grida, gridb, gridc)
            self.kpoints["shift"] = shift
        else:
            with open(kpts, "r") as infile:
                kp_list = infile.readlines().split()
            self.kpoints["numkp"] = numkp
            self.kpoints["points"] = np.array(kp_list, dtype='d')

        return

    # TODO: atomic position change and adding/deleting atoms based on the space group operators
    def unitcell_transform(self, unit=None, unitvec=None, rotate=None):
        if unitvec is not None:
            self.cellparam["vector"] = np.array(unitvec, ntype='d')
            self.cellparam["unit"] = unit

        if rotate is not None:
            self.cellparam["vector"] = np.dot(self.cellparam["vector"], np.array(rotate, dtype='d'))

        if unit != self.cellparam["unit"]:
            newvec = []
            for x in self.cellparam["vector"]:
                newvec.append(Unitconverter.unit_convert(x, "length", self.cellparam["unit"], unit))
            self.cellparam["vector"] = np.array(newvec, dtype='d')
            self.cellparam["unit"] = unit

        return

    # TODO: handling alat and crystal_sg units
    def position_unittransform(self, unit):
        if unit != self.atompos["unit"]:
            newpos = []
            if unit == "crystal":
                for x in self.atompos["coordinates"]:
                    newpos.append(np.dot(x, np.linalg.inv(self.cellparam["vector"])))
            else:
                for x in self.atompos["coordinates"]:
                    newpos.append(Unitconverter.unit_convert(x, "length", self.cellparam["unit"], unit))
            self.atompos["coordinates"] = np.array(newpos, dtype='d')
        else:
            pass

        return

    def position_translation(self, transvec, transunit):
        newpos = []
        if transunit != self.atompos["unit"]:
            if self.atompos["unit"] == "crystal":
                transvec = np.dot(transvec, np.linalg.inv(self.cellparam["vector"]))
            else:
                transvec = Unitconverter.unit_convert(transvec, "length", transunit, self.sellparam["unit"])
        else:
            transvec = transvec

        for x in self.atompos["coordinates"]:
            newpos.append(x + np.array(transvec, dtype='d'))

        self.atompos["coordinates"] = np.array(newpos, dtype='d')
        return

    def add_atoms(self, elements, coordinates, unit):
        if len(elements) != len(coordinates):
            raise IOError("Number of input elements and coordinates not matching!")

        for x in elements:
            self.atompos["elements"] = np.append(self.atompos["elements"], x)

        if unit != self.atompos["unit"]:
            newcoord = []
            if self.atompos["unit"] == "crystal":
                for x in coordinates:
                    newcoord.append(np.dot(x, np.linalg.inv(self.cellparam["vector"])))
            else:
                for x in coordinates:
                    newpos.append(Unitconverter.unit_convert(x, "length", unit, self.sellparam["unit"]))
        else:
            pass

        for x in np.array(newcoord, dtype='d'):
            self.atompos["coordinates"] = np.append(self.atompos["coordinates"], x)

        self.atompos["coordinates"] = np.reshape(self.atompos["coordinates"], (len(self.atompos["elements"]), 3))
        return


def tag_parser(linein):
    line = linein.partition("=")
    tag = str(line[0]).strip()
    val = str(line[2]).strip()
    return tag, val


def coord_parser(linein):
    elem = linein.split()[0]
    coord = np.array([linein.split()[1], linein.split()[2], linein.split()[3]], dtype='d')
    if len(linein.split()) > 4:
        dynamics = np.array([linein.split()[4], linein.split()[5], linein.split()[6]], dtype='d')
    else:
        dynamics = np.array([None, None, None])
    return elem, coord, dynamics


def pseudo_parser(linein):
    dic = {}
    elem = linein.split()[0]
    mass = linein.split()[1]
    pp = linein.split()[2]
    dic["mass"] = mass
    dic["pseudo"] = pp
    return elem, dic


def cell_volume(unitvec):
    volume = np.dot(np.cross(unitvec[1], unitvec[2]), unitvec[0])
    return volume


def rec_vector(unitvec):
    a_star = np.cross(unitvec[1], unitvec[2]) / cell_volume(unitvec)
    b_star = np.cross(unitvec[2], unitvec[0]) / cell_volume(unitvec)
    c_star = np.cross(unitvec[0], unitvec[1]) / cell_volume(unitvec)
    return np.array([a_star, b_star, c_star], dtype='d')


def ordered_load(stream, yamlloader=yaml.Loader, object_pairs_hook=OrderedDict):
    class OrderedLoader(yamlloader):
        pass

    # noinspection PyArgumentList
    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))

    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)
    return yaml.load(stream, OrderedLoader)
