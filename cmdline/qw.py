#!/usr/local/bin/python

import argparse
import sys
import numpy as np
from pwscf.pw_input import PWin
from qwutils.plotigor import PlotIgor


def executepwintags(args):
    p = PWin()
    p.read_file(args.input)

    if args.write is not None:
        p.write_tags_from_tagvaluepairlist(args.write)

    if args.remove is not None:
        p.remove_tags(args.remove)

    p.write_pwin(args.output)
    return


def executepwincell(args):
    p = PWin()
    p.read_file(args.input)
    p.unitcell_transform(args.unit, args.cell, args.rotation)
    p.write_pwin(args.output)
    return


def executepwinatom(args):
    p = PWin()
    p.read_file(args.input)

    if args.elements is not None:
        p.write_atoms(args.elements, args.coordinates, args.unit)

    else:
        if args.unit is not None:
            if args.translation is not None:
                p.position_translation(args.translation, args.unit)
            else:
                p.position_unittransform(args.unit)

    p.write_pwin(args.output)
    return


def executepwinpseudo(args):
    p = PWin()
    p.read_file(args.input)

    if args.pseudomass is None:
        mass = np.full(len(args.pseudoelem), None)
    else:
        mass = args.pseudomass

    if args.pseudofile is not None:
        for i in range(len(args.pseudoelem)):
            p.write_pseudo(args.pseudoelem[i], mass[i], args.pseudofile[i])
    else:
        for i in range(len(args.pseudoelem)):
            p.write_pseudo(args.pseudoelem[i], mass[i], None, args.pseudolabel[i])

    p.write_pwin(args.output)
    return


def executepwinkpoints(args):
    p = PWin()
    p.read_file(args.input)
    p.change_kpoints(args.ktype, args.kspacing, args.kgrid, args.kshift, args.knum, args.kpts, args.kreduce)
    p.write_pwin(args.output)
    return


def executepwingen(args):
    p = PWin()
    p.read_file(args.input)
    p.generate_pwin(args.tags)
    p.write_pwin(args.output)
    return


def executeplotband(args):
    p = PlotIgor(args.input, args.output, args.prefix)
    print("Reading %s ... " % args.input)
    p.read_bands()
    print("Done!")
    print("Writing %s ... " % args.output)
    p.write_bands(args.plot, args.fermi, args.shift, args.guide)
    print("Done!")
    return


def executeplotdos(args):
    p = PlotIgor(args.input, args.output, args.prefix)
    print("Reading %s ... " % args.input)
    p.read_dos()
    print("Done!")
    print("Writing %s ... " % args.output)
    p.write_dos(args.plot, args.fermi)
    print("Done!")
    return


def executeplotwf(args):
    p = PlotIgor(args.input, args.output, args.prefix)
    print("Reading %s ... " % args.input)
    p.read_wf()
    print("Done!")
    print("Writing %s ... " % args.output)
    p.write_wf(args.plot)
    print("Done!")
    return


def executeplotdiel(args):
    p = PlotIgor(None, args.output, args.prefix)
    print("Parsing dielectric functions ... ")
    p.read_diel(args.real, args.imag, args.diag, args.eels, args.direction)
    print("Done!")
    print("Writing %s ... " % args.output)
    p.write_diel(args.plot)
    print("Done!")
    return


def main():
    description = """qw.py
    Type qw.py [function] --help for more detailed informations.
    """

    desc_in = """
    """

    # Main parser
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

    # Subparsers
    subparsers = parser.add_subparsers(title="Functions")

    parser_pwin = subparsers.add_parser("pwin", formatter_class=argparse.RawTextHelpFormatter, description=desc_in)

    pwinsubparsers = parser_pwin.add_subparsers()

    parser_tags = pwinsubparsers.add_parser("tags")
    parser_tags.add_argument("-i", dest="input", type=str, default="pw.in")
    parser_tags.add_argument("-o", dest="output", type=str, default="pw.in")
    parser_tags.add_argument("-w", dest="write", type=str, default=None, nargs='*')
    parser_tags.add_argument("-r", dest="remove", type=str, default=None, nargs='*')
    parser_tags.set_defaults(func=executepwintags)

    parser_cell = pwinsubparsers.add_parser("cell")
    parser_cell.add_argument("-i", dest="input", type=str, default="pw.in")
    parser_cell.add_argument("-o", dest="output", type=str, default="pw.in")
    parser_cell.add_argument("-u", dest="unit", type=str, default=None)
    parser_cell.add_argument("-v", dest="cell", type=str, default=None, nargs='*')
    parser_cell.add_argument("-r", dest="rotation", type=str, default=None, nargs='*')
    parser_cell.set_defaults(func=executepwincell)

    parser_atom = pwinsubparsers.add_parser("atom")
    parser_atom.add_argument("-i", dest="input", type=str, default="pw.in")
    parser_atom.add_argument("-o", dest="output", type=str, default="pw.in")
    parser_atom.add_argument("-e", dest="elements", type=str, default=None, nargs='*')
    parser_atom.add_argument("-c", dest="coordinates", type=str, default=None, nargs='*')
    parser_atom.add_argument("-u", dest="unit", type=str, default=None)
    parser_atom.add_argument("-t", dest="translation", type=str, default=None, nargs='*')
    parser_atom.set_defaults(func=executepwinatom)

    parser_pseudo = pwinsubparsers.add_parser("pseudo")
    parser_pseudo.add_argument("-i", dest="input", type=str, default="pw.in")
    parser_pseudo.add_argument("-o", dest="output", type=str, default="pw.in")
    parser_pseudo.add_argument("-pe", dest="pseudoelem", type=str, default=None, nargs='*')
    parser_pseudo.add_argument("-pm", dest="pseudomass", type=str, default=None, nargs='*')
    parser_pseudo.add_argument("-pl", dest="pseudolabel", type=str, default=None, nargs='*')
    parser_pseudo.add_argument("-pf", dest="pseudofile", type=str, default=None, nargs='*')
    parser_pseudo.set_defaults(func=executepwinpseudo)

    parser_kpoints = pwinsubparsers.add_parser("kpoints")
    parser_kpoints.add_argument("-ktype", dest="ktype", type=str, default=None)
    parser_kpoints.add_argument("-kgrid", dest="kgrid", type=str, default=None, nargs='*')
    parser_kpoints.add_argument("-kshift", dest="kshift", type=str, default=None, nargs='*')
    parser_kpoints.add_argument("-kspacing", dest="kspacing", type=str, default=None)
    parser_kpoints.add_argument("-knum", dest="knum", type=str, default=None)
    parser_kpoints.add_argument("-klist", dest="klist", type=str, default=None, nargs='*')
    parser_kpoints.add_argument("-kred", dest="kreduce", type=str, default=None)
    parser_kpoints.set_defaults(func=executepwinkpoints)

    parser_gen = pwinsubparsers.add_parser("generate")
    parser_gen.add_argument("-t", dest="tags", type=str, default=None, nargs='*')
    parser_gen.set_defaults(func=executepwingen)

    parser_plot = subparsers.add_parser("plot", formatter_class=argparse.RawTextHelpFormatter, description=desc_in)

    plotsubparsers = parser_plot.add_subparsers()

    parser_band = plotsubparsers.add_parser("band")
    parser_band.add_argument("-i", dest="input", type=str, default="bands.dat")
    parser_band.add_argument("-o", dest="output", type=str, default="bands.itx")
    parser_band.add_argument("-p", dest="prefix", type=str, default=None)
    parser_band.add_argument("-f", dest="fermi", type=float, default=0.0)
    parser_band.add_argument("-s", dest="shift", action='store_true')
    parser_band.add_argument("-g", dest="guide", action='store_true')
    parser_band.add_argument("-P", dest="plot", action='store_false')
    parser_band.set_defaults(func=executeplotband)

    parser_dos = plotsubparsers.add_parser("dos")
    parser_dos.add_argument("-i", dest="input", type=str, default="bands.dat")
    parser_dos.add_argument("-o", dest="output", type=str, default="bands.itx")
    parser_dos.add_argument("-p", dest="prefix", type=str, default=None)
    parser_dos.add_argument("-f", dest="fermi", type=float, default=0.0)
    parser_dos.add_argument("-P", dest="plot", action='store_false')
    parser_dos.set_defaults(func=executeplotdos)

    parser_pdos = plotsubparsers.add_parser("pdos")
    parser_pdos.add_argument("-i", dest="input", type=str, default="pdos_tot")
    parser_pdos.add_argument("-o", dest="output", type=str, default="pdos.itx")
    parser_pdos.add_argument("-p", dest="prefix", type=str, default=None)
    parser_pdos.add_argument("-f", dest="fermi", type=float, default=0.0)
    parser_pdos.add_argument("-P", dest="plot", action='store_false')
    parser_pdos.add_argument("-c", dest="combine", action='store_false')
    parser_pdos.set_defaults(func=executeplotpdos)

    parser_wf = plotsubparsers.add_parser("wf")
    parser_wf.add_argument("-i", dest="input", type=str, default="avg.dat")
    parser_wf.add_argument("-o", dest="output", type=str, default="avg.itx")
    parser_wf.add_argument("-p", dest="prefix", type=str, default=None)
    parser_wf.add_argument("-P", dest="plot", action='store_false')
    parser_wf.set_defaults(func=executeplotwf)

    parser_diel = plotsubparsers.add_parser("diel")
    parser_diel.add_argument("-i", dest="imag", type=str, default="epsi_pwscf.dat")
    parser_diel.add_argument("-r", dest="real", type=str, default="epsr_pwscf.dat")
    parser_diel.add_argument("-d", dest="diag", type=str, default="ieps_pwscf.dat")
    parser_diel.add_argument("-e", dest="eels", type=str, default="eels_pwscf.dat")
    parser_diel.add_argument("-o", dest="output", type=str, default="epsilon.itx")
    parser_diel.add_argument("-p", dest="prefix", type=str, default=None)
    parser_diel.add_argument("-P", dest="plot", action='store_false')
    parser_diel.add_argument("-D", dest="direction", action='store_true')
    parser_diel.set_defaults(func=executeplotdiel)

    args = parser.parse_args()

    try:
        getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)


if __name__ == "__main__":
    main()
