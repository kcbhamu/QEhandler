from __future__ import unicode_literals, print_function
from collections import OrderedDict
from qwutils.generalutils import ScientificConstants, Unitconverter
import numpy as np
import yaml
import spglib


class PWout(object):

    def __init__(self, infile, calculation=None):
        self.infile = infile
        self.memory = None
        self.calculation = None or calculation
        self.status = None
        self.energy = None
        self.structure = None
        return

    def read_file(self):
        with open(self.infile, "r") as infile:
            f = infile.readlines()
            # for line in f:
            #     return

        return
