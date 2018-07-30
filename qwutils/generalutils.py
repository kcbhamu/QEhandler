from ase import data, units


class Unitconverter(object):
    def __init__(self):
        return

    @staticmethod
    def unit_convert(value, system, unitin, unitout):
        dic = {"energy": {"Ry": units.Ry, "Rydberg": units.Ry, "rydberg": units.Ry, "eV": units.eV,
                          "electronvolt": units.eV, "J": units.J, "Joule": units.J, "joule": units.J,
                          "Ha": units.Ha, "Hartree": units.Ha, "hartree": units.Ha},
               "length": {"Ang": units.Ang, "angstrom": units.Ang, "ang": units.Ang, "nanometer": units.nm,
                          "nm": units.nm, "Nanometer": units.nm, "m": units.m, "meter": units.m, "Meter": units.m,
                          "Bohr": units.Bohr, "bohr": units.Bohr},
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
