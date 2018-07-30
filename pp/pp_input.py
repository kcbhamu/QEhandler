from __future__ import unicode_literals, print_function
from collections import OrderedDict
from qwutils.generalutils import ScientificConstants, Unitconverter
import os
import sys
import yaml
import numpy as np

class PPin(object):

    def __init__(self):
        self.inputpp = OrderedDict()
        self.plot = OrderedDict()
        with open(os.path.join(os.path.dirname(__file__), "pp_tags.yaml"), "r") as y:
            self.taglist = ordered_load(y, yaml.SafeLoader)[0]

