import pdb
import sys
import json
import os.path
import glob, os
import csv

# additional package
from itertools import permutations
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
import pandas as pd
import numpy as np


data_dir = './core/data'



# add paths for novostoic
sys.path.append('./core/')
#sys.path.append('./core/data')


#from reprime import *
#from novoStoic import *




def all_perms(elements):
    if len(elements) <=1:
        yield elements
    else:
        for perm in all_perms(elements[1:]):
            for i in range(len(elements)):
                # nb elements[0:1] works in both string and list contexts
                yield perm[:i] + elements[0:1] + perm[i:]


def intermediate_tracer(substrates,products,novel_mets,solutions,\
                radius,reaction_dict,primary_substrate,primary_product,project):
    
    met_save = dict()
    met_result = {}
    cofactors = pd.read_csv('./core/data/cofactors.csv',index_col=0)


    molsigs = pd.read_csv(os.path.join(data_dir,'molecular_signature_' + str(radius) + '.csv'),\
            index_col=0).fillna(0)
    if not novel_mets:
        molsigs_product = molsigs[primary_product]
    else:
        molsigs_product = pd.read_csv('./test/novel_mets/molsig_'+ primary_product[0] + '_' + \
                    str(radius) + '.csv', index_col=0)
    novel_moiet = list(molsigs_product.index)
    rules = pd.read_csv(os.path.join(data_dir,'rule_' + str(radius) + '.csv'),index_col=0)
    moiety_index = rules.index.tolist()

    reaction_list = list(pd.read_csv('./core//data/rule_0.csv', index_col=0).columns.values)
    exchange_index = substrates.keys() + products.keys()
    T = rules.to_dict(orient='index')
    C = molsigs.to_dict(orient='index')
    #molsigs_product = molsigs_product.to_dict(orient='index')

    for i in moiety_index:
        if i in novel_moiet:
            check = 1
        else:
            #pdb.set_trace()
            df = pd.DataFrame([float(0)],columns = primary_product)
            s = df.xs(0)
            s.name = i
            molsigs_product = molsigs_product.append(s)


    C_novel = molsigs_product.to_dict(orient='index')
    sol_iter = 0
    sol_save = dict()
    #pdb.set_trace()
    store_met = 0
    for things in solutions:
        #db.set_trace()
        sol_num = 0
        met_save[store_met] = dict()
        sol_save[store_met] = dict()
        sol_rules = solutions[things]
    #sol_rules = {'R05283': -2,'R00713': -1,'R01173': -1}
    
        sol_names = sol_rules.keys()
        #['R05283','R00713','R01173']

        T2 = dict()

        for m in moiety_index:
            T2[m] = dict()
        for sol in sol_rules:

            for m in moiety_index:
                T2[m][sol] = T[m][sol]
            for metab in reaction_dict[sol]:
                #pdb.set_trace()
                if metab in cofactors.iloc[:,0]:
                    for m in moiety_index:
                        T2[m][sol] = T2[m][sol] - reaction_dict[sol][metab]*C[m][metab]
 
 
        for sol in sol_rules:
            if abs(sol_rules[sol]) > 1:
                for i in range(1,int(abs(sol_rules[sol]))+1):
                    for m in moiety_index:
                        T2[m][sol + '_' + str(int(abs(sol_rules[sol])))] = T2[m][sol]
        for sol in sol_rules.keys():
            if abs(sol_rules[sol]) > 1:
                for i in range(2,int(abs(sol_rules[sol]))+1):
                    sol_rules[sol + '_' + str(i)] = sol_rules[sol]
                    sol_names.append(sol + '_' + str(i))
        perm = all_perms(sol_names)
        sol_perm = list(perm)
        #pdb.set_trace()
        minim = dict()
        
        store = 0
        #pdb.set_trace()
        for i in range(0,len(sol_perm)):
            sol_rules_check = sol_perm[i]
            sol = sol_rules_check[0]
            Rupdate = dict()
            met = dict()
        #next step: incorporate permutations into sorting 

            met[0] = dict()
            for m in moiety_index:
                #pdb.set_trace()
                met[0][m] = C[m][primary_substrate[0]]
            met[1] = dict()
            k = 1
            #pdb.set_trace()
            proceed = 1
            for m in moiety_index:
                if proceed == 1:
                    if sol_rules[sol] < 0:
                        Rupdate[m] = T2[m][sol] + substrates[primary_substrate[0]] * C[m][primary_substrate[0]]
                        Rupdate[m] = -Rupdate[m]
                    elif sol_rules[sol] > 0:
                        Rupdate[m] = T2[m][sol] - substrates[primary_substrate[0]] * C[m][primary_substrate[0]]
                        #Rupdate[m] = -Rupdate[m]
                    met[1][m] = Rupdate[m]
                    if met[1][m] < 0:
                        #pdb.set_trace()
                        proceed = 0
            if proceed == 1:
                for sol2 in sol_rules_check:
                    if sol2 == sol:
                        for m in moiety_index:
                            Rupdate[m] = Rupdate[m] 
                    else:
                        k = k + 1
                        met[k] = dict()
                        for m in moiety_index:
                            if proceed == 1:
                                if sol_rules[sol2] < 0:
                                    Rupdate[m] = T2[m][sol2] - Rupdate[m]
                                    Rupdate[m] = -Rupdate[m]
                                    #pdb.set_trace()
                                elif sol_rules[sol2] >0:
                                    Rupdate[m] = T2[m][sol2] + Rupdate[m]
                                    #Rupdate[m] = -Rupdate[m]
                                    #pdb.set_trace()
                            met[k][m] = Rupdate[m]
                            if met[k][m] < 0:
                                proceed = 0
                if proceed == 1:
                    met_save[store_met][store] = dict()
                    k = k + 1
                    met[k] = dict()
                    for m in moiety_index:
                        Rupdate[m] = Rupdate[m] - C_novel[m][primary_product[0]]
                        met[k][m] = C_novel[m][primary_product[0]]                
                    storage_k = 0
                    for ind in range( store * k , store * k + k ):
                    #for ind in range( 0, len(met.keys())):
                        met_save[store_met][store][storage_k] = dict()
                        for m in moiety_index:
                            met_save[store_met][store][storage_k][m] = met[ind - store*k][m]
                            #met_save[store_met][store][storage_k][m] = met[ind][m]
                        storage_k = storage_k + 1
                    for key,value in met_save.items():
                        if value not in met_result.values():
                            met_result[key] = value
                    
                    store = store + 1        
                
                    
                    #Rcheck = 0
                    #pdb.set_trace()
                    #for m in moiety_index:
                    #    Rcheck = abs(Rupdate[m]) + Rcheck    
                    #if Rcheck == 0:
                    #with open(project + '/result/'+'/sol_order_'+str(sol_iter)+'_'+str(sol_num)+'.csv','wb') as f:
                    #    f.write('Reaction order: \n')
                    #    for j in sol_rules_check:
                    #        f.write(j)
                    #        f.write('\n')
                    
                    sol_save[store_met][sol_num] = sol_rules_check
                    sol_num = sol_num + 1
                        #pdb.set_trace()
        store_met = store_met + 1
        sol_iter = sol_iter + 1
        
    #for i in met_save.keys():
    #    if not met_save[i]:
    #        del met_save[i]

            #pdb.set_trace()
        #with open('intermediate_metabolites.json', 'w') as fp:
        #    json.dump(met_result, fp)

        #check_semi = dict()
#        for i in met_result:
#            check_semi[i] = dict()
#            check_semi[i] = 0
#            for j in met_result[i]:
#                check_semi[i] = check_semi[i] + abs(met_result[i][j] - molsigs[primary_substrate][j])
    remove = dict()

    for i in sol_save:
        for k in sol_save[i]:
            for j in range(0,len(sol_save[i][k])):
                sep = '_'
                bb = sol_save[i][k][j]
                bb2 = bb.split(sep,1)[0]
                sol_save[i][k][j] = bb2
    #pdb.set_trace()
    for i in sol_save:
        remove[i] = dict()
        #pdb.set_trace()
        for k in sol_save[i]:
            sol_enum = 0
            for jj in sol_save[i]:
                check_iden = cmp(sol_save[i][k], sol_save[i][jj])
                if check_iden == 0:
                    sol_enum = sol_enum + 1
                    if sol_enum > 1 and jj == k:
                        remove[i][k] = max(jj,i)
                #pdb.set_trace()
    for i in remove:
        #pdb.set_trace()
        remove[i] = sorted(remove[i], key=remove[i].get, reverse=True)
        for j in range(0,len(remove[i])):
            del sol_save[i][remove[i][j]]
            del met_result[i][remove[i][j]]

                
    #pdb.set_trace()
    return (met_save, sol_save)

