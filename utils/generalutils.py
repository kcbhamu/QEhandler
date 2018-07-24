from ase import data, units


class Unitconverter(object):
    def __init__(self):
        return

    @staticmethod
    def unit_convert(value, system, unitin, unitout):
        dic = {"energy": {"Ry": units.Ry, "eV": units.eV, "J": units.J, "Ha": units.Ha},
               "length": {"Ang": units.Ang, "nm": units.nm, "m": units.m, "Bohr": units.Bohr},
               }
        try:
            return value * (dic[system][unitin] / dic[system][unitout])
        except KeyError:
            raise IOError("Error: Unknown unit.")


class ScientificConstants(object):
    def __init__(self):
        return

    @staticmethod
    def elementary_data(element, key=None):
        number = data.atomic_numbers[element]

        if key is None:
            value = data.atomic_masses[number]
        elif key == "mass":
            value = data.atomic_masses[number]
        elif key == "vdw":
            value = data.vdw_radii[number]
        elif key == "covalent":
            value = data.covalent_radii[number]

        return value
