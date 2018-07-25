#!/usr/local/bin/python

import argparse
import os
import sys
import numpy as np
from pwscf.pw_input import PWin

# installpath = '/Users/Woosun/Dropbox/Dev/QEhandler'
# sys.path.extend([installpath])


def executepwintags(args):
    p = PWin(args.input, args.output)
    p.read_file()

    if args.write is not None:
        p.write_tags_from_tagvaluepairlist(args.write)

    if args.remove is not None:
        p.remove_tags(args.remove)

    p.write_pwin(args.output)
    return


def executepwinstructure(args):
    p = PWin(args.input, args.output)
    p.read_file()

    return


def executepwinpseudo(args):
    p = PWin(args.input, args.output)
    p.read_file()

    return


def executepwinkpoints(args):
    p = PWin(args.input, args.output)
    p.read_file()

    return


def executepwingen(args):
    p = PWin(args.input, args.output)
    p.generate_pwin()
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

    parser_structure = pwinsubparsers.add_parser("structure")
    parser_structure.add_argument("-i", dest="input", type=str, default="pw.in")
    parser_structure.add_argument("-o", dest="output", type=str, default="pw.in")
    parser_structure.add_argument("-cell", dest="cell", type=str, default=None, nargs='*')
    parser_structure.add_argument("-e", dest="elements", type=str, default=None, nargs='*')
    parser_structure.add_argument("-c", dest="coordinates", type=str, default=None, nargs='*')
    parser_structure.add_argument("-cu", dest="coordinateunit", type=str, default=None, nargs='*')
    parser_structure.add_argument("-r", dest="rotation", type=str, default=None, nargs='*')
    parser_structure.add_argument("-t", dest="translation", type=str, default=None, nargs='*')
    parser_structure.set_defaults(func=executepwinstructure)

    parser_pseudo = pwinsubparsers.add_parser("pseudo")
    parser_pseudo.add_argument("-i", dest="input", type=str, default="pw.in")
    parser_pseudo.add_argument("-o", dest="output", type=str, default="pw.in")
    parser_pseudo.add_argument("-pe", dest="pseudoelem", type=str, default=None, nargs='*')
    parser_pseudo.add_argument("-pm", dest="pseudomass", type=str, default=None, nargs='*')
    parser_pseudo.add_argument("-pl", dest="pseudolabel", type=str, default=None, nargs='*')
    parser_pseudo.add_argument("-pf", dest="pseudofile", type=str, default=None, nargs='*')
    parser_pseudo.set_defaults(func=executepwinpseudo)

    parser_kpoints = pwinsubparsers.add_parser("kpoints")
    parser_kpoints.add_argument("-ktype", dest="ktype", type=str, default=None, nargs='*')
    parser_kpoints.add_argument("-kgrid", dest="kgrid", type=str, default=None, nargs='*')
    parser_kpoints.add_argument("-kshift", dest="kshift", type=str, default=None, nargs='*')
    parser_kpoints.add_argument("-kspacing", dest="kspacing", type=str, default=None, nargs='*')
    parser_kpoints.add_argument("-knum", dest="knum", type=str, default=None, nargs='*')
    parser_kpoints.add_argument("-klist", dest="klist", type=str, default=None, nargs='*')
    parser_kpoints.add_argument("-kred", dest="kreduce", type=str, default=None, nargs='*')
    parser_kpoints.set_defaults(func=executepwinkpoints)

    parser_gen = pwinsubparsers.add_parser("generate")
    parser_gen.add_argument("-i", dest="input", type=str, default="pw.in")
    parser_gen.add_argument("-o", dest="output", type=str, default="pw.in")
    parser_gen.add_argument("-d", dest="default", action='store_true')
    parser_gen.add_argument("-t", dest="tags", type=str, default=None, nargs='*')
    parser_gen.add_argument("-cell", dest="cell", type=str, default=None, nargs='*')
    parser_gen.add_argument("-e", dest="elements", type=str, default=None, nargs='*')
    parser_gen.add_argument("-c", dest="coordinates", type=str, default=None, nargs='*')
    parser_gen.add_argument("-cu", dest="coordinateunit", type=str, default=None, nargs='*')
    parser_gen.add_argument("-pe", dest="pseudoelem", type=str, default=None, nargs='*')
    parser_gen.add_argument("-pm", dest="pseudomass", type=str, default=None, nargs='*')
    parser_gen.add_argument("-pl", dest="pseudolabel", type=str, default=None, nargs='*')
    parser_gen.add_argument("-pf", dest="pseudofile", type=str, default=None, nargs='*')
    parser_gen.add_argument("-kgrid", dest="kgrid", type=str, default=None, nargs='*')
    parser_gen.add_argument("-kshift", dest="kshift", type=str, default=None, nargs='*')
    parser_gen.add_argument("-kspacing", dest="kspacing", type=str, default=None, nargs='*')
    parser_gen.add_argument("-knum", dest="knum", type=str, default=None, nargs='*')
    parser_gen.add_argument("-klist", dest="klist", type=str, default=None, nargs='*')
    parser_gen.add_argument("-kred", dest="kreduce", type=str, default=None, nargs='*')
    parser_gen.set_defaults(func=executepwingen)

    args = parser.parse_args()

    try:
        getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)


if __name__ == "__main__":
    main()
