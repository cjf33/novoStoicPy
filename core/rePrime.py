import numpy
import cobra
import pdb
import sys
import json
import re
import os.path
import glob, os
sys.path.append('../')

from util import *
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops, AddHs
import pandas as pd

data_dir = './core/data'
test_dir = './test/novel_mets'

def byteify(input):
    if isinstance(input, dict):
        return {byteify(key): byteify(value)
                for key, value in input.iteritems()}
    elif isinstance(input, list):
        return [byteify(element) for element in input]
    elif isinstance(input, unicode):
        return input.encode('utf-8')
    else:
        return input

def count_mol(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

# calculate molecular signature
def init_with_atom_feature(metList,init_val):
    molecular_signature = dict()
    all_features = dict()
    all_features[0] = set()
    i = 0
    missing = []
    for index, met in metList.iterrows():
        #if met in missing: continue
        #if met == 'C00080': continue # this is a proton, ignore it
        #if met == 'C00282': continue # this is h2, ignore it
        #pdb.set_trace()
        try:
            # molecular_signature = dict()
            mol = Chem.MolFromSmiles(met['SMILES'])
            G = mol_to_nx(mol)
            #print(G)
            # assign feature string                             
            met_features, G = init_atomfeature(G)
            all_features[0] = all_features[0].union(met_features)
            #print(all_features)
            #pdb.set_trace()
            nx.write_yaml(G, os.path.join(data_dir,'radius_0/' + met['MNX_ID'] + '.yaml'))
        except Exception as e:
            print met['MNX_ID']
            missing.append(met['MNX_ID'])
            i = i + 1
            #pdb.set_trace()
   #pdb.set_trace()
    if init_val == 1:
        df = pd.DataFrame(sorted(all_features[0]),columns=None)
        df.to_csv('moieties_0.csv',header=None,index=False)
        mdf = pd.DataFrame(missing,columns=None)
        mdf.to_csv('missing.csv',header=None,index=False)
    else:
        molsigs = pd.read_csv(os.path.join(data_dir,'moieties_0.csv'),header=None)
        df = pd.DataFrame(sorted(all_features[0]),columns=None)
        for i in range(0,len(df)):
            df_check = molsigs.isin([df.iloc[i,0]])
            if not any(df_check.iloc[:,0]):
                molsigs = molsigs.append([df.iloc[i,0]],ignore_index=True)
                molsigs.to_csv(os.path.join(data_dir,'moieties_0.csv'),header=None,index=False)
        #pdb.set_trace()
def get_prime(metList,radius):
    for index, met in metList.iterrows():
        #pdb.set_trace()
        if pd.isnull(met['InChI']): continue
        else:
            #pdb.set_trace()
        #if met in missing: continue
        #if met == 'C00080': continue # this is just a proton, ignore it
        #if met == 'C00282': continue # this is just h2, ignore it
            try:      
            # at radius = 0, moieties
            # radius = 0
                G = nx.read_yaml(os.path.join(data_dir,'./radius_' +str(radius) + '/' + met['MNX_ID'] + '.yaml')) 
                product_list = pd.read_csv(os.path.join(data_dir,'./moieties_' + str(radius) + '.csv'),\
                                                header=None)[0].tolist()
                #pdb.set_trace()
                if radius >= 1:
                    product_list = map(int, product_list)
                G = processPrm(G,radius,product_list)
                nx.write_yaml(G,os.path.join(data_dir,'./radius_' +str(radius) + '/' + met['MNX_ID'] + '.yaml'))
            except Exception as e:
                print met['MNX_ID']
                raise e

def get_product(metList,radius,init_val):
    molecular_signature = dict()
    all_features = dict()
    all_features[0] = set()
    for index, met in metList.iterrows():
        #if met in missing: continue
        #if met == 'C00080': continue # this is just a proton, ignore it
        #if met == 'C00282': continue # this is just h2, ignore it
        if pd.isnull(met['InChI']): continue
        else:
            try:
                # pdb.set_trace()
                G = nx.read_yaml(os.path.join(data_dir,'./radius_' +str(radius-1) + '/' + met['MNX_ID'] + '.yaml')) 

                # assign feature string                             
                met_features, G = processPrd(G,radius)
                all_features[0] = all_features[0].union(met_features)
                # pdb.set_trace()
                nx.write_yaml(G,os.path.join(data_dir,'./radius_' +str(radius) + '/' + met['MNX_ID'] + '.yaml'))
            except Exception as e:
                print met
                #pdb.set_trace()
                raise e
    if init_val == 1:
        df = pd.DataFrame(sorted(all_features[0]),columns=None)
        df.to_csv(os.path.join(data_dir,'moieties_'+str(radius)+'.csv'),header=None,index=False)
        mdf = pd.DataFrame(missing,columns=None)
        mdf.to_csv(os.path.join(data_dir,'missing.csv'),header=None,index=False)
    else:
        molsigs = pd.read_csv(os.path.join(data_dir,'moieties_'+str(radius)+'.csv'),header=None)
        df = pd.DataFrame(sorted(all_features[0]),columns=None)
        for i in range(0,len(df)):
            df_check = molsigs.isin([df.iloc[i,0]])
            if not any(df_check.iloc[:,0]):
                molsigs = molsigs.append([df.iloc[i,0]],ignore_index=True)
                molsigs.to_csv(os.path.join(data_dir,'moieties_'+str(radius)+'.csv'),header=None,index=False)
        #pdb.set_trace()

# get molecular signature: count the number of moieties
def calmolecularsignature(metList,radius):
    molecular_signature = dict()
    for index, met in metList.iterrows():
        if pd.isnull(met['InChI']): continue
        else:
            #if met in missing: continue
            #if met == 'C00080': continue # this is just a proton, ignore it
            #if met == 'C00282': continue # this is just h2, ignore it
            G = nx.read_yaml('./radius_' +str(radius) + '/' + met['MNX_ID'] + '.yaml')
            # reactantsGraph.append(primesGraph)
            cardinality_dict = cal_cardinality(G,radius)
            # reactantsCardinality.append(cardinality_dict)
            molecular_signature[met['MNX_ID']] = cardinality_dict
    molsigna_df = pd.DataFrame.from_dict(molecular_signature).fillna(0)
    molsigna_df.to_csv('molecular_signature_update_'+str(radius) + '.csv',index=True)
    pdb.set_trace()
    return molecular_signature
def calunknownsignature(metList):
    #moieties = pd.read_csv(os.path.join(data_dir,'moieties_'+str(radius)+'.csv'),header=None)
    #primes = loadPrimes()
    #prime_update = [None]*(len(moieties))
    #mol = MolFromSmiles()
    for index, met in metList.iterrows():
        #G = nx.read_yaml(os.path.join(test_dir,met['MNX_ID']+'.yaml'))
        mol = Chem.MolFromSmiles(met['SMILES'])
        #mol = AddHs(mol)
        G = mol_to_nx(mol)
        # assign feature string                             
        met_features,G = init_atomfeature(G)
        molecular_signature = dict()
        #cardinality_dict = cal_cardinality(G,0)
        #molecular_signature_ = pd.DataFrame.from_dict(cardinality_dict,orient='index',columns=[met['MNX_ID']])
        #molecular_signature = molecular_signature.fillna(0)
        #molecular_signature.to_csv(os.path.join(test_dir,'molsig_'+met+'_'+str(radius) + '.csv'),index=True)

        for radius in range(0,2):
            product_list = pd.read_csv('./core/data/moieties_' + str(radius) + '.csv',\
                                                    header=None)[0].tolist()
            processPrm(G,radius,product_list)
            met_features, G = processPrd(G,radius+1)
        #nx.write_yaml(G,filename)
        for radius in range(0,3):
            cardinality_dict = cal_cardinality(G,radius)
            molecular_signature[met['ID']] = cardinality_dict
            molsigna_df = pd.DataFrame.from_dict(molecular_signature).fillna(0)
            molsigna_df.to_csv('./test/novel_mets/molsig___' + met['ID'] + '_' + str(radius) +\
                             '.csv',index=True)



def updatemolecularsignature(metList,radius):
    molecular_signature = pd.read_csv(os.path.join(data_dir,'molecular_signature_'+ str(radius)+'.csv'),index_col=0)
    rules1 = pd.read_csv(os.path.join(data_dir,'rule_'+ str(radius)+'_noduplic.csv'),index_col=0)
    rules2 = pd.read_csv(os.path.join(data_dir,'rule_'+ str(radius)+'.csv'),index_col=0)
    moieties = pd.read_csv(os.path.join(data_dir,'moieties_'+str(radius)+'.csv'),header=None)
    primes = loadPrimes()
    prime_update = [None]*(len(moieties)-len(molecular_signature))
    if len(moieties) != len(molecular_signature):

        for i in range(0,len(moieties)-len(molecular_signature)):
            prime_update[i] = primes[len(molecular_signature) + i]
        #pdb.set_trace()
        update_mol_s = numpy.zeros([len(moieties)-len(molecular_signature),len(molecular_signature.columns)])
        update_rules1 = numpy.zeros([len(moieties)-len(molecular_signature),len(rules1.columns)])
        update_rules2 = numpy.zeros([len(moieties)-len(molecular_signature),len(rules2.columns)])
        molecular_signature = molecular_signature.append(pd.DataFrame(update_mol_s,columns=molecular_signature.columns.tolist(),index=prime_update))
        rules1 = rules1.append(pd.DataFrame(update_rules1,columns=rules1.columns.tolist(),index=prime_update))
        rules2 = rules2.append(pd.DataFrame(update_rules2,columns=rules2.columns.tolist(),index=prime_update))

    for index, met in metList.iterrows():
        if pd.isnull(met['InChI']): continue
        else:
            G = nx.read_yaml(os.path.join(data_dir,'./radius_' +str(radius) + '/' + met['MNX_ID'] + '.yaml'))
            cardinality_dict = cal_cardinality(G,radius)
            molecular_signature_new = pd.DataFrame.from_dict(cardinality_dict,orient='index',columns=[met['MNX_ID']])
            molecular_signature = molecular_signature.add(molecular_signature_new,fill_value=0)
    molecular_signature = molecular_signature.fillna(0)
    molecular_signature.to_csv(os.path.join(data_dir,'molecular_signature_'+str(radius) + '.csv'),index=True)
    rules1.to_csv(os.path.join(data_dir,'rule_'+str(radius) + '_noduplic.csv'),index=True)
    rules2.to_csv(os.path.join(data_dir,'rule_'+str(radius) + '.csv'),index=True)
    return molecular_signature

def update_rxn_rule(newrule,radius):
    mol_sig = pd.read_csv('./core/data/molecular_signature_'+str(radius)+'.csv',index_col=0)

    rules2 = pd.read_csv(os.path.join(data_dir,'rule_'+ str(radius)+'.csv'),index_col=0)
    sig_vals = mol_sig.index.values.tolist()
    mol_sig = mol_sig.to_dict()
    rxn_list = newrule.to_dict()['Formula']
    rxn_update = dict()
    jsonFile = open("./core/data/optstoic_v3_Sji_dict.json", "r+")
    data = json.load(jsonFile)

    if radius == 0:
        for mnx in rxn_list.keys():
            if mnx not in data.keys():
                with open("./core/data/optstoic_v3_reduced_reactions.txt","a") as f:
                    f.write("\n"+mnx)

    for mnx in rxn_list.keys():
        if mnx not in data.keys():
            rxn_update[mnx] = dict()
            rl_stoic = rxn_list[mnx].split(" + ")
            matching = [s for s in rl_stoic if "<=>" in s]
            matching2 = rl_stoic.index(matching[0])
            for i in range(len(rl_stoic)):
                if i < matching2:
                    r_pull = rl_stoic[i].split(")")
                    r_pull1 = float(r_pull[0].replace("(",""))
                    r_pull2 = r_pull[1].replace(" ","")
                    rxn_update[mnx][r_pull2] = -r_pull1
                elif i == matching2:
                    r_pull = rl_stoic[i].split(" <=> ")
                    r_pull1 = r_pull[0].split(")")
                    r_pull11 = float(r_pull1[0].replace("(",""))
                    r_pull12 = r_pull1[1].replace(" ","")
                    rxn_update[mnx][r_pull12] = -r_pull11
                    r_pull2 = r_pull[1].split(")")
                    r_pull21 = float(r_pull2[0].replace("(",""))
                    r_pull22 = r_pull2[1].replace(" ","")
                    rxn_update[mnx][r_pull22] = r_pull11
                else:
                    r_pull = rl_stoic[i].split(")")
                    r_pull1 = float(r_pull[0].replace("(",""))
                    r_pull2 = r_pull[1].replace(" ","")
                    rxn_update[mnx][r_pull2] = r_pull1
            data[mnx] = rxn_update[mnx]
    with open("./core/data/optstoic_v3_Sji_dict.json", "w") as outfile:
        json.dump(data,outfile)
    rule_new = dict()
    rule_empty = dict()
    for rxn in rxn_update:
        rule_new[rxn] = dict()
        for moiety in sig_vals:
            rule_new[rxn][moiety] = 0 
        for metab in rxn_update[rxn].keys():
            for moiety in sig_vals:
                rule_new[rxn][moiety] = rule_new[rxn][moiety] + rxn_update[rxn][metab]*mol_sig[metab][moiety]

    rules_check = pd.DataFrame.from_dict(rule_new,orient='columns')
    rules_update2 = pd.concat([rules2,rules_check],axis=1)
    rules_update1 = rules_update2.T.drop_duplicates().T
    rules_update1.to_csv(os.path.join(data_dir,'rule_'+str(radius) + '_noduplic.csv'),index=True)
    rules_update2.to_csv(os.path.join(data_dir,'rule_'+str(radius) + '.csv'),index=True)
    #pdb.set_trace()

# calculate reaction rule
def get_rxn_rule(radius):
    reaction_dict = json.load(open('./test.json','rb'))
    molsigna_df = pd.read_csv('molecular_signature_'+ str(radius)+'.csv',index_col=0)
    missing = pd.read_csv('missing.csv',names=['metab'])
    #missing = missing.metab.tolist()
    #pdb.set_trace()
    rule_df = pd.DataFrame(index=molsigna_df.index)
    for rid,value in reaction_dict.iteritems():
        value=byteify(value)
        # skip the reactions with missing metabolites
        #b = value.keys()
        #pdb.set_trace()
        if not set(missing['metab']).isdisjoint(value.keys()): continue

        rule_df[rid] = 0
        for met,stoic in value.iteritems():
            if met == 'MNXM01' or met == 'MNXM1': continue # hydogen is zero
            rule_df[rid] += molsigna_df[met]*stoic

    rule_df.to_csv('rule_'+ str(radius)+'.csv',index=True)
    #pdb.set_trace()
def remove_duplicate(radius):
    rule1 = pd.read_csv('rule_'+str(radius) + '.csv',index_col=0).T.drop_duplicates().T
    rule1.to_csv('rule_'+str(radius) + '_noduplic.csv',index=True)

def get_molsig_exmetab(metList):
    
    for index, met in metList.iterrows():
        mol = Chem.MolFromSmiles(met['SMILES'])
        #mol = Chem.MolFromMolFile('../../novoStoic_data/KEGG/' + met + '.mol')
        filename = './test/novel_mets/' + met['ID'] + '.yaml'

        G = mol_to_nx(mol)
        # assign feature string                             
        met_features,G = init_atomfeature(G)
        molecular_signature = dict()

        for radius in range(0,3):
            product_list = pd.read_csv('./core/data/moieties_' + str(radius) + '.csv',\
                                                    header=None)[0].tolist()
            processPrm(G,radius,product_list)
            met_features, G = processPrd(G,radius+1)
        nx.write_yaml(G,filename)

        for radius in range(0,3):
            cardinality_dict = cal_cardinality(G,radius)
            molecular_signature[met['ID']] = cardinality_dict
            molsigna_df = pd.DataFrame.from_dict(molecular_signature).fillna(0)
            molsigna_df.to_csv('./test/novel_mets/molsig_'+ met['ID'] + '_' + str(radius) +\
                                 '.csv',index=True)


