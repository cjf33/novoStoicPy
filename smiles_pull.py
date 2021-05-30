import sys
import pdb
import os.path
import os
import json
import pandas as pd
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import AddHs

import PIL
from PIL import Image
from PIL import ImageChops
import numpy as np
import graphviz as gv
sys.path.append('./core/')
from structure_gen import *
data_dir = './core/data'
radius = 1



rules = {0: ['MNXR92214', 'MNXR86251', 'MNXR92214'], 1: ['MNXR86251', 'MNXR92214', 'MNXR92214']}
solution = {'MNXR92214': -2.0, 'MNXR92214_2': -2.0, 'MNXR86251': -1.0}
metab_dict = {0: {0: u'MNXM92', 1: 'MNXM514', 2: 'MNXM2334', 3: u'ChEBI_41189'}, 1: {0: u'MNXM92', 1: 'Hypothetical_Metabolite_0', 2: 'MNXM2334', 3: u'ChEBI_41189'}}
r_type_input = {0: {0: 'rule', 1: 'rule', 2: 'rule'}, 1: {0: 'rule', 1: 'rule', 2: 'rule'}}
iteration = 0
r_type = {'MNXR92214': 'rule', 'MNXR86251': 'rule'}
project = u'test'
product_name = u'14bdo'
out_ind = 0
reaction_dict = json.load(open('./core/data/metanetx_sji.json','rb'))
substrate = [u'MNXM92']


for i in rules:
	(out_ind) = structure_gen(rules,substrate,solution,iteration,r_type,metab_dict,project,reaction_dict,r_type_input,product_name,out_ind)
	out_ind = out_ind + 1

