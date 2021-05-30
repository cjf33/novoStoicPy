import pdb
import sys
import json
import os.path
import glob, os

# additional package
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
import pandas as pd

# add paths for novostoic
sys.path.append('./core/')
#sys.path.append('./core/data')
from novoStoic2 import *

name = 'bdo'

##### novoStoic
novoStoic(name)