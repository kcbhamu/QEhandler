from __future__ import unicode_literals, print_function
from collections import OrderedDict
import numpy as np
import yaml


# TODO: writing all namelist tags at the end of class calling to prevent overlapping writing methods
class PWin(object):

    def __init__(self, file, cellparam=None, atompos=None, aspecies=None, kpoints=None, control=None, system=None,
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
        self.infile = file or {}

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

            if "CELL_PARAMETERS" in namelist:
                dic["CELL_PARAMETERS"]["vector"] = np.array(cell_param, dtype='d')
                self.cellparam = dic["CELL_PARAMETERS"]

            if "K_POINTS" in namelist:
                if dic["K_POINTS"]["unit"] == "automatic":
                    dic["K_POINTS"]["grid"] = [kp_list[0][0], kp_list[0][1], kp_list[0][2]]
                    dic["K_POINTS"]["shift"] = [kp_list[0][3], kp_list[0][4], kp_list[0][5]]
                elif dic["K_POINTS"]["unit"] == "gamma":
                    pass
                else:
                    dic["K_POINTS"]["numpoints"] = kp_list.pop(0)
                    dic["K_POINTS"]["points"] = np.array(kp_list, dtype='d')
                self.kpoints = dic["K_POINTS"]

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

    def generate_pwin(self, tagvaluepair=None):
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


def ordered_load(stream, yamlloader=yaml.Loader, object_pairs_hook=OrderedDict):
    class OrderedLoader(yamlloader):
        pass

    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))

    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)
    return yaml.load(stream, OrderedLoader)

