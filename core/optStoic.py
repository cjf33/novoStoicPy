#/usr/bin/python

# author: Lin Wang
"""
 optStoic: Designing overall stoichiometric conversions and intervening 
 metabolic reactions

 https://doi.org/10.1038/srep16009

"""
import pulp
import json
import re
import pandas as pd
import os
import pdb
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops


#pulp_solver = pulp.solvers.CPLEX_CMD(path=None, keepFiles=0, mip=1, msg=1,\
#      options=['mip tolerances mipgap 0', 'mip tolerances absmipgap 0',\
#       'mip tolerances integrality 0', 'simplex tolerances optimality 1E-9',\
#       'simplex tolerances feasibility 1E-9',], timelimit=1200)

pulp_solver = pulp.solvers.CPLEX_PY()
#GUROBI_CMD_OPTIONS = [('Threads', 2), ('TimeLimit', 1200), \
#                ('MIPGapAbs', 1e-6), ('MIPGap', 1e-6), ('CliqueCuts', 2)]
#pulp_solver = pulp.solvers.GUROBI_CMD(path=None, keepFiles=0, mip=1, msg=1,
#                options=GUROBI_CMD_OPTIONS)
# pulp_solver = pulp.solvers.GLPK_CMD(path=None, keepFiles=0, mip=1, msg=1, options=[])
data_dir = './core/data'

def construct_metS(list_of_mets,name):
    all_mets_S = json.load(open(os.path.join(data_dir, "met_details_dict_v2.json")))
    #pdb.set_trace()
    elements = ['C','H','N','O','P','S','F','Cl','Mg','Fe','Se','Co','As','Br',
                'I','R','charge']
    atomic_number = [6, 1, 7, 8, 15, 16, 9, 17, 12, 26, 34, 27, 33, 35, 53, 500000]

    #pdb.set_trace()
    missingMet = pd.read_excel(name+"_optstoic_input.xlsx",sheet_name='missingMet',header=None,names=['MNX_ID','InChI','dGf'])
    miss_met = dict()
    for index1, met in missingMet.iterrows():

        if met['MNX_ID'] not in all_mets_S:
            mol = Chem.MolFromInchi(met['InChI'].encode('ascii', 'ignore'))
            mol = Chem.AddHs(mol)
            #pdb.set_trace()
            miss_met[met['MNX_ID']] = dict()
            for index in range(0,len(atomic_number)):
                num_atoms = count_mol(mol,atomic_number[index])
                miss_met[met['MNX_ID']][elements[index]] = float(num_atoms)
                #pdb.set_trace()
            miss_met[met['MNX_ID']]['charge'] = float(Chem.GetFormalCharge(mol))
            miss_met[met['MNX_ID']]['dGf'] = float(met['dGf'])
            #miss_met[met['MNX_ID']]['H'] = Chem.AddHs(mol) - mol.GetNumAtoms()
            #pdb.set_trace()
    all_mets_S.update(miss_met)
    #pdb.set_trace()


    met_S = {met_id: all_mets_S[met_id] for met_id in list_of_mets}
    #pdb.set_trace()
    return met_S


def count_mol(mol,atomic_number):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == atomic_number)

def simulate_optStoic(name):
    dG = pd.read_csv(name+"_dG_output.csv")
    reactant_in = pd.read_excel(open(name+'_optstoic_input.xlsx','rb'),sheet_name='reactant_stoichs',header=None,index_col=0)
    reactant_in = reactant_in[1].to_dict()
    for i in reactant_in.keys():
        reactant_in[i] = None if reactant_in[i] == 'None' else reactant_in[i]
    product_in = pd.read_excel(open(name+'_optstoic_input.xlsx','rb'),sheet_name='product_stoichs',header=None,index_col=0)
    product_in = product_in[1].to_dict()
    for i in product_in.keys():
        product_in[i] = None if product_in[i] == 'None' else product_in[i]
    objective_in = pd.read_excel(open(name+'_optstoic_input.xlsx','rb'),sheet_name='objective',header=None,index_col=0)
    objective_in = objective_in.index.tolist()[0]
    params = {}
    params['reactant_stoichs'] = reactant_in
    params['product_stoichs'] = product_in
    params['objective'] = objective_in
    list_of_mets = []
    print params['reactant_stoichs']
    for met in params['reactant_stoichs']:
        list_of_mets.append(met)
    for met in params['product_stoichs']:
        list_of_mets.append(met)

    print list_of_mets

    # add proton to the metabolit list because many time it is required
    if 'MNXM1' not in list_of_mets:
        list_of_mets.append('MNXM1')
    met_S = construct_metS(list_of_mets,name)
    for i in range(0,len(dG)):
        met_S[dG.loc[i,'Compund']]['dGf'] = dG.loc[i,'dG (kJ/mol)'] 
    metabolites = met_S.keys()
    dGmax = 5
    M = 100
    # objective = 'MaxTargetYield'
    EPS = 1e-5
    # allow_list_metabolite = ['C00007',#    /*o2*/
    #                         'C00267',#    /*glc*/
    #                         # 'C00011',#    /*co2*/
    #                         'C00080',#    /*h+*/
    #                         'C00001',#    /*h2o*/
    #                         # 'C02457',#    /*13pdo*/
    #                         'C00033',#    /*acetate*/
    #                         ]

    #------- define variables
    # zp      objective function
    # s(i)    stoichiometric coefficient of metabolite i
    # ;
    s = pulp.LpVariable.dicts("s", metabolites, lowBound=-M, upBound=M, \
                                cat='Continuous')


    #------- define LP problem
    lp_prob = pulp.LpProblem("OptStoic", pulp.LpMaximize)
    

    #------- define objective function
    # obj..                     zp =e= sum(i$pdt(i),s(i))/(-s('%substrate%'));
    objective_stoic = params['objective']
    lp_prob += s[objective_stoic], "MaxTargetYield"
    #pdb.set_trace()

    #------- define constraints
    # stoic(j)$(elem(j))..   sum(i,s(i)*m(i,j)) =e= 0;
    # bal(j)$(dG(j))..      sum(i,s(i)*m(i,j)) =l= dGmax*(-s('%substrate%'));

    elements = ['C','H','N','O','P','S','F','Cl','Mg','Fe','Se','Co','As','Br',
                'I','R','charge']
    #atom_num = [6, 1, 7, 8, 15, 16, 9, 17, 12, 26, 34, 27, 33, 35, 53]


    for j in elements:
        lp_prob += pulp.lpSum(s[i]*met_S[i][j] for i in metabolites) == 0

    lp_prob += pulp.lpSum(s[i]*met_S[i]['dGf'] for i in metabolites) <= \
                dGmax,'element_'+ j
                # dGmax*s[substrate_metabolite],'element_'+ j

    for met,stoic in params['reactant_stoichs'].iteritems():
        # if data['fixed_stoich'] != None:
        #     stoic = -data['fixed_stoich']
        #     lp_prob += s[data['start_compound_id'][0]] == stoic
        if stoic != None:
            lp_prob += s[met] == stoic
        else:
            lp_prob += s[met] <= 0
    for met,stoic in params['product_stoichs'].iteritems():
        # if data['fixed_stoich2'] != None:
        #     stoic = data['fixed_stoich2']
        #     lp_prob += s[data['target_compound_id'][0]] == stoic
        if stoic != None:
            lp_prob += s[met] == stoic
        else:
            lp_prob += s[met] >= 0


    #------- solve the problem
    lp_prob.solve(pulp_solver)

    Ex_stoic = {}
    for i in metabolites:
        if s[i].varValue is not None:
            if s[i].varValue > EPS or s[i].varValue < -EPS:
                print i, s[i].varValue
                Ex_stoic[i] = s[i].varValue

    stoic_out = pd.DataFrame.from_dict(Ex_stoic, orient="index")
    stoic_out.to_csv(name+"_optstoic_output.csv",header=None)
    return Ex_stoic

if __name__ == '__main__':
    params = {}

    # Example 1: glucose -> acetate 
    #params['reactant_stoichs'] = {
    #    'C00267': -1, # glc
    #
    #    }
    #params['product_stoichs'] = {
    #    'C00033': None, # acetate
    #    }
    #
    #params['objective'] = 'C00033'
    #Ex_stoic = simulate_optStoic(params)
    #print Ex_stoic

    # exampel 2:
    params['reactant_stoichs'] = {
        'C00091': -1, # succo
        'C00005': -4, # NADPH
        'C00080': None, # H+
        }
    params['product_stoichs'] = {
        'ChEBI_41189': 1, # 1,4-butanediol
        'C00006': 4, # NADP+
        'C00010': 1, # CoA
        'C00001': 1, # H20
        }
    params['objective'] = 'ChEBI_41189'
    #.set_trace()
    Ex_stoic = simulate_optStoic(params)
    print Ex_stoic   