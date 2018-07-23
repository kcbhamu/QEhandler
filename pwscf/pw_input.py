from __future__ import unicode_literals, print_function
from collections import OrderedDict
import numpy as np
import yaml


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

    def generate_draft(self, tags, values):
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

        cell = [("", "")]

        self.control = OrderedDict(control)
        self.system = OrderedDict(system)
        self.electrons = OrderedDict(electrons)
        self.cell = OrderedDict(cell)

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

    def remove_tags(self, tag, value):
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


def ordered_load(stream, Loader=yaml.Loader, object_pairs_hook=OrderedDict):
    class OrderedLoader(Loader):
        pass
    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))
    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)
    return yaml.load(stream, OrderedLoader)


# taglist = OrderedDict({"control": {"calculation": {"values": ["scf", "nscf", "bands", "relax",
#                                                               "md", "vc-relax", "vc-md"],
#                                                    "default": "scf",
#                                                    "type": "str"},
#                                    "title": {"values": None,
#                                              "default": "",
#                                              "type": "str"},
#                                    "verbosity": {"values": ["low", "high"],
#                                                  "default": "low",
#                                                  "type": "str"},
#                                    "restart_mode": {"values": ["from_scratch", "restart"],
#                                                     "default": "from_scratch",
#                                                     "type": "str"},
#                                    "wf_collect": {"default": ".TRUE.",
#                                                   "type": "logical"},
#                                    "nstep": {"default": "1",
#                                              "default2": "50",
#                                              "type": "int"},
#                                    "iprint": {"default": None,
#                                               "type": "int"},
#                                    "tstress": {"default": ".false.",
#                                                "type": "logical"},
#                                    "tprnfor": {"default": None,
#                                                "type": "logical"},
#                                    "dt": {"default": "20.D0",
#                                           "type": "real"},
#                                    "outdir": {"default": "./",
#                                               "type": "str"},
#                                    "wfcdir": {"default": "./",
#                                               "type": "str"},
#                                    "prefix": {"default": "pwscf",
#                                               "type": "str"},
#                                    "lkpoint_dir": {"default": ".true.",
#                                                    "type": "logical"},
#                                    "max_seconds": {"default": "1.D+7",
#                                                    "type": "real"},
#                                    "etot_conv_thr": {"default": "1.0D-4",
#                                                      "type": "real"},
#                                    "forc_conv_thr": {"default": "1.0D-3",
#                                                      "type": "real"},
#                                    "disk_io": {"values": ["high", "medium", "low", "none"],
#                                                "default": "low",
#                                                "default2": "medium",
#                                                "type": "str"},
#                                    "pseudo_dir": {"default": None,
#                                                   "type": "str"},
#                                    "tefield": {"default": ".FALSE.",
#                                                "type": "logical"},
#                                    "dipfield": {"default": ".FALSE.",
#                                                 "type": "logical"},
#                                    "lelfield": {"default": ".FALSE.",
#                                                 "type": "logical"},
#                                    "nberrycyc": {"default": "1",
#                                                  "type": "int"},
#                                    "lorbm": {"default": ".FALSE.",
#                                              "type": "logical"},
#                                    "lberry": {"default": ".FALSE.",
#                                               "type": "logical"},
#                                    "gdir": {"values": ["1", "2", "3"],
#                                             "default": None,
#                                             "type": "int"},
#                                    "nppstr": {"default": None,
#                                               "type": "int"},
#                                    "lfcpopt": {"default": ".FALSE.",
#                                                "type": "logical"},
#                                    "gate": {"default": ".FALSE.",
#                                             "type": "logical"},
#                                    },
#                        "system": {"ibrav": {"values": ["0", "1", "2", "3", "-3", "4", "5",
#                                                        "-5", "6", "7", "8", "9", "-9", "10",
#                                                        "11", "12", "-12", "13", "14"],
#                                             "default": "0",
#                                             "type": "int",
#                                             "required": True
#                                             },
#                                   "celldm": {"default": None,
#                                              "type": "real"
#                                              },
#                                   "A": {"default": None,
#                                         "type": "real"
#                                         },
#                                   "B": {"default": None,
#                                         "type": "real"
#                                         },
#                                   "C": {"default": None,
#                                         "type": "real"
#                                         },
#                                   "cosAB": {"default": None,
#                                             "type": "real"
#                                             },
#                                   "cosAC": {"default": None,
#                                             "type": "real"
#                                             },
#                                   "cosBC": {"default": None,
#                                             "type": "real"
#                                             },
#                                   "nat": {"default": None,
#                                           "type": "int",
#                                           "required": True
#                                           },
#                                   "ntyp": {"default": None,
#                                            "type": "int",
#                                            "required": True
#                                            },
#                                   "nbnd": {"default": None,
#                                            "type": "int"
#                                            },
#                                   "tot_charge": {"default": "0.0",
#                                                  "type": "real"
#                                                  },
#                                   "starting_charge": {"default": "0.0",
#                                                       "type": "real"
#                                                       },
#                                   "tot_magnetization": {"default": "-1",
#                                                         "type": "real"
#                                                         },
#                                   "starting_magnetization(i)": {"default": None,
#                                                                 "type": "real"
#                                                                 },
#                                   "ecutwfc": {"default": None,
#                                               "type": "real",
#                                               "required": True
#                                               },
#                                   "ecutrho": {"default": None,
#                                               "type": "real",
#                                               },
#                                   "ecutfock": {"default": None,
#                                                "type": "real",
#                                                },
#                                   "nr1": {"default": None,
#                                           "type": "int",
#                                           },
#                                   "nr2": {"default": None,
#                                           "type": "int",
#                                           },
#                                   "nr3": {"default": None,
#                                           "type": "int",
#                                           },
#                                   "nr1s": {"default": None,
#                                            "type": "int",
#                                            },
#                                   "nr2s": {"default": None,
#                                            "type": "int",
#                                            },
#                                   "nr3s": {"default": None,
#                                            "type": "int",
#                                            },
#                                   "nosym": {"default": ".FALSE.",
#                                             "type": "logical",
#                                             },
#                                   "nosym_evc": {"default": ".FALSE.",
#                                                 "type": "logical",
#                                                 },
#                                   "noinv": {"default": ".FALSE.",
#                                             "type": "logical",
#                                             },
#                                   "no_t_rev": {"default": ".FALSE.",
#                                                "type": "logical",
#                                                },
#                                   "force_symmorphic": {"default": ".FALSE.",
#                                                        "type": "logical",
#                                                        },
#                                   "use_all_frac": {"default": ".FALSE.",
#                                                    "type": "logical",
#                                                    },
#                                   "occupations": {
#                                       "values": ["smearing", "tetrahedra", "tetrahedra_lin", "tetrahedra_opt",
#                                                  "fixed", "from_input"],
#                                       "default": None,
#                                       "type": "str",
#                                       },
#                                   "one_atom_occupations": {"default": ".FALSE.",
#                                                            "type": "logical",
#                                                            },
#                                   "starting_spin_angle": {"default": ".FALSE.",
#                                                           "type": "logical",
#                                                           },
#                                   "degauss": {"default": "0.D0",
#                                               "type": "real",
#                                               },
#                                   "smearing": {"values": ["gaussian", "gauss", "methfessel-paxton", "m-p", "mp",
#                                                           "mazari-vanderbilt", "cold", "m-v", "mv",
#                                                           "fermi-dirac", "f-d", "fd"],
#                                                "default": "gaussian",
#                                                "type": "str",
#                                                },
#                                   "nspin": {"default": "1",
#                                             "type": "int",
#                                             },
#                                   "noncolin": {"default": ".FALSE.",
#                                                "type": "logical",
#                                                },
#                                   "ecfixed": {"default": "0.0",
#                                               "type": "real",
#                                               },
#                                   "qcutz": {"default": "0.0",
#                                             "type": "real",
#                                             },
#                                   "q2sigma": {"default": "0.1",
#                                               "type": "real",
#                                               },
#                                   "input_dft": {"values": [],
#                                                 "default": None,
#                                                 "type": "str",
#                                                 },
#                                   "exx_fraction": {"default": None,
#                                                    "type": "real",
#                                                    },
#                                   "screening_parameter": {"default": "0.106",
#                                                           "type": "real",
#                                                           },
#                                   "exxdiv_treatment": {
#                                       "values": ["gygi-baldereschi", "vcut_spherical", "vcut_ws", "none"],
#                                       "default": "gygi-baldereschi",
#                                       "type": "str",
#                                       },
#                                   "x_gamma_extrapolation": {"default": ".true.",
#                                                             "type": "logical",
#                                                             },
#                                   "ecutvcut": {"default": "0.0",
#                                                "type": "real",
#                                                },
#                                   "nqx1": {"default": None,
#                                            "type": "int",
#                                            },
#                                   "nqx2": {"default": None,
#                                            "type": "int",
#                                            },
#                                   "nqx3": {"default": None,
#                                            "type": "int",
#                                            },
#                                   "lda_plus_u": {"default": ".FALSE.",
#                                                  "type": "logical",
#                                                  },
#                                   "lda_plus_u_kind": {"default": "0",
#                                                       "type": "int",
#                                                       },
#                                   "Hubbard_U": {"default": None,
#                                                 "type": "real",
#                                                 },
#                                   "Hubbard_J0": {"default": None,
#                                                  "type": "real",
#                                                  },
#                                   "Hubbard_alpha": {"default": None,
#                                                     "type": "real",
#                                                     },
#                                   "Hubbard_beta": {"default": None,
#                                                    "type": "real",
#                                                    },
#                                   "Hubbard_J(i,ityp)": {"default": None,
#                                                         "type": "real",
#                                                         },
#                                   "starting_ns_eigenvalue(m,ispin,l)": {"default": "-1.D0",
#                                                                         "type": "real",
#                                                                         },
#                                   "U_projection_type": {
#                                       "values": ["atomic", "ortho-atomic", "norm-atomic", "file", "pseudo"],
#                                       "default": "atomic",
#                                       "type": "str",
#                                       },
#                                   "edir": {"default": None,
#                                            "type": "int",
#                                            },
#                                   "emaxpos": {"default": "0.5D0",
#                                               "type": "real",
#                                               },
#                                   "eopreg": {"default": "0.5D0",
#                                              "type": "real",
#                                              },
#                                   "eamp": {"default": "0.001",
#                                            "type": "real",
#                                            },
#                                   "angle1": {"default": None,
#                                              "type": "real",
#                                              },
#                                   "angle2": {"default": None,
#                                              "type": "real",
#                                              },
#                                   "constrained_magnetization": {"values": ["none", "total", "atomic",
#                                                                            "total direction", "atomic direction"],
#                                                                 "default": "none",
#                                                                 "type": "str",
#                                                                 },
#                                   "fixed_magnetization": {"default": "0.D0",
#                                                           "type": "real",
#                                                           },
#                                   "lambda": {"default": "1.D0",
#                                              "type": "real",
#                                              },
#                                   "report": {"default": "100",
#                                              "type": "int",
#                                              },
#                                   "lspinorb": {"default": None,
#                                                "type": "logical",
#                                                },
#                                   "assume_isolated": {"values": ["none", "makov-payne", "m-p", "mp",
#                                                                  "martyna-tuckerman", "m-t", "mt", "esm"],
#                                                       "default": "none",
#                                                       "type": "str",
#                                                       },
#                                   "esm_bc": {"values": ["pbc", "bc1", "bc2", "bc3"],
#                                              "default": "pbc",
#                                              "type": "str",
#                                              },
#                                   "esm_w": {"default": "0.D0",
#                                             "type": "real",
#                                             },
#                                   "esm_efield": {"default": "0.D0",
#                                                  "type": "real",
#                                                  },
#                                   "esm_nfit": {"default": "4",
#                                                "type": "int",
#                                                },
#                                   "fcp_mu": {"default": "0.D0",
#                                              "type": "real",
#                                              },
#                                   "vdw_corr": {"values": ["grimme-d2", "Grimme-D2", "DFT-D", "dft-d",
#                                                           "TS", "ts", "ts-vdw", "ts-vdW", "tkatchenko-scheffler",
#                                                           "XDM", "xdm"],
#                                                "default": "none",
#                                                "type": "str",
#                                                },
#                                   "london": {"default": ".FALSE.",
#                                              "type": "logical",
#                                              },
#                                   "london_s6": {"default": "0.75",
#                                                 "type": "real",
#                                                 },
#                                   "london_c6": {"default": None,
#                                                 "type": "real",
#                                                 },
#                                   "london_rvdw": {"default": None,
#                                                   "type": "real",
#                                                   },
#                                   "london_rcut": {"default": "200",
#                                                   "type": "real",
#                                                   },
#                                   "ts_vdw_econv_thr": {"default": "1.D-6",
#                                                        "type": "real",
#                                                        },
#                                   "ts_vdw_isolated": {"default": ".FALSE.",
#                                                       "type": "logical",
#                                                       },
#                                   "xdm": {"default": ".FALSE.",
#                                           "type": "logical",
#                                           },
#                                   "xdm_a1": {"default": "0.6836",
#                                              "type": "real",
#                                              },
#                                   "xdm_a2": {"default": "1.5045",
#                                              "type": "real",
#                                              },
#                                   "space_group": {"default": "0",
#                                                   "type": "int",
#                                                   },
#                                   "uniqueb": {"default": ".FALSE.",
#                                               "type": "logical",
#                                               },
#                                   "origin_choice": {"default": "1",
#                                                     "type": "int",
#                                                     },
#                                   "rhombohedral": {"default": ".TRUE.",
#                                                    "type": "logical",
#                                                    },
#                                   "zgate": {"default": "0.5",
#                                             "type": "real",
#                                             },
#                                   "relaxz": {"default": ".FALSE.",
#                                              "type": "logical",
#                                              },
#                                   "block": {"default": ".FALSE.",
#                                             "type": "logical",
#                                             },
#                                   "block_1": {"default": "0.45",
#                                               "type": "real",
#                                               },
#                                   "block_2": {"default": "0.55",
#                                               "type": "real",
#                                               },
#                                   "block_height": {"default": "0.1",
#                                                    "type": "real",
#                                                    },
#
#                                   },
#                        "electrons": {"electron_maxstep": {"default": "100",
#                                                           "type": "int",
#                                                           },
#                                      "scf_must_converge": {"default": ".TRUE.",
#                                                            "type": "logical",
#                                                            },
#                                      "conv_thr": {"default": "1.D-6",
#                                                   "type": "real",
#                                                   },
#                                      "adaptive_thr": {"default": ".FALSE.",
#                                                       "type": "logical",
#                                                       },
#                                      "conv_thr_init": {"default": "1.D-3",
#                                                        "type": "real",
#                                                        },
#                                      "conv_thr_multi": {"default": "1.D-1",
#                                                         "type": "real",
#                                                         },
#                                      "mixing_mode": {"values": ["plain", "TF", "local-TF"],
#                                                      "default": "plain",
#                                                      "type": "str",
#                                                      },
#                                      "mixing_beta": {"default": "0.7D0",
#                                                      "type": "real",
#                                                      },
#                                      "mixing_ndim": {"default": "8",
#                                                      "type": "int",
#                                                      },
#                                      "mixing_fixed_ns": {"default": "0",
#                                                          "type": "int",
#                                                          },
#                                      "diagonalization": {"values": ["david", "cg", "cg-serial", "david-serial"],
#                                                          "default": "david",
#                                                          "type": "str",
#                                                          },
#                                      "ortho_para": {"default": "0",
#                                                     "type": "int",
#                                                     },
#                                      "diago_thr_init": {"default": None,
#                                                         "type": "real",
#                                                         },
#                                      "diago_cg_maxiter": {"default": None,
#                                                           "type": "int",
#                                                           },
#                                      "diag_david_ndim": {"default": "4",
#                                                          "type": "int",
#                                                          },
#                                      "diago_full_acc": {"default": ".FALSE.",
#                                                         "type": "logical",
#                                                         },
#                                      "efield": {"default": "0.D0",
#                                                 "type": "real",
#                                                 },
#                                      "efield_cart": {"default": "(0.D0, 0.D0, 0.D0)",
#                                                      "type": "real",
#                                                      },
#                                      "efield_phase": {"values": ["read", "write", "none"],
#                                                       "default": "none",
#                                                       "type": "str",
#                                                       },
#                                      "startingpot": {"values": ["atomic", "file"],
#                                                      "default": None,
#                                                      "type": "str",
#                                                      },
#                                      "startingwfc": {"values": ["atomic", "atomic+random", "random", "file"],
#                                                      "default": None,
#                                                      "type": "str",
#                                                      },
#                                      "tqr": {"default": ".FALSE.",
#                                              "type": "logical",
#                                              },
#                                      },
#                        "ions": {
#                            "ion_dynamics": {"values": ["bfgs", "damp", "verlet", "langevin", "langevin-smc", "beeman"],
#                                             "default": None,
#                                             "type": "str",
#                                             },
#                            "ion_positions": {"values": ["default", "from_input"],
#                                              "default": "default",
#                                              "type": "str",
#                                              },
#                            "pot_extrapolation": {"values": ["none", "atomic", "first_order", "second_order"],
#                                                  "default": "atomic",
#                                                  "type": "str",
#                                                  },
#                            "wfc_extrapolation": {"values": ["none", "first_order", "second_order"],
#                                                  "default": "none",
#                                                  "type": "str",
#                                                  },
#                            "remove_rigid_rot": {"default": ".FALSE.",
#                                                 "type": "logical",
#                                                 },
#                            "ion_temperature": {"values": ["rescaling", "rescale-v", "rescale-T", "reduce-T",
#                                                           "berendsen", "andersen", "initial", "not_controlled"],
#                                                "default": "not_controlled",
#                                                "type": "str",
#                                                },
#                            "tempw": {"default": "300.D0",
#                                      "type": "real",
#                                      },
#                            "tolp": {"default": "100.D0",
#                                     "type": "real",
#                                     },
#                            "delta_t": {"default": "1.D0",
#                                        "type": "real",
#                                        },
#                            "nraise": {"default": "1",
#                                       "type": "int",
#                                       },
#                            "refold_pos": {"default": ".FALSE.",
#                                           "type": "logical",
#                                           },
#                            "upscale": {"default": "100.D0",
#                                        "type": "real",
#                                        },
#                            "bfgs_ndim": {"default": "1",
#                                          "type": "int",
#                                          },
#                            "trust_radius_max": {"default": "0.8D0",
#                                                 "type": "real",
#                                                 },
#                            "trust_radius_min": {"default": "1.D-3",
#                                                 "type": "real",
#                                                 },
#                            "trust_radius_ini": {"default": "0.5D0",
#                                                 "type": "real",
#                                                 },
#                            "w_1": {"default": "0.01D0",
#                                    "type": "real",
#                                    },
#                            "w_2": {"default": "0.5D0",
#                                    "type": "real",
#                                    },
#                            },
#                        "cell": {"cell_dynamics": {"values": ["none", "sd", "damp-pr", "damp-w", "bfgs", "pr", "w"],
#                                                   "default": None,
#                                                   "type": "str",
#                                                   },
#                                 "press": {"default": "0.D0",
#                                           "type": "real",
#                                           },
#                                 "wmass": {"default": None,
#                                           "type": "real",
#                                           },
#                                 "cell_factor": {"default": "1.0",
#                                                 "default2": "2.0",
#                                                 "type": "real",
#                                                 },
#                                 "press_conv_thr": {"default": "0.5D0",
#                                                    "type": "real",
#                                                    },
#                                 "cell_dofree": {"values": ["all", "x", "y", "z", "xy", "xz", "yz", "xyz",
#                                                            "shape", "volume", "2Dxy", "2Dshape"],
#                                                 "default": "all",
#                                                 "type": "str",
#                                                 }
#                                 }
#                        })
