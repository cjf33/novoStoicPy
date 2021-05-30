import pandas as pd
import pulp
import pdb
import itertools
import json
import os
from rePrime import *
from intermediate_tracer import *
from structure_gen import *

missing = []

data_dir = './core/data'


def get_S(reaction_dict):
    """build the stoichiometric matrix at a specific growth rate"""
    # reaction_dict = json.load(open('./optstoic_v3_reduced_Sji.json','rb'))
    S = {}

    # populate with stoichiometry
    for r, stoic_info in reaction_dict.iteritems():
        # if r == "R00402": pdb.set_trace()
        b = stoic_info.keys()
        if not set(missing).isdisjoint(b): continue

        for met, value in stoic_info.iteritems():
            if met not in S:
                S[met] = {}
                S[met][r] = float(value)
            else:
                S[met][r] = float(value)
    #pdb.set_trace()
    return S


def novoStoic(name):
    reaction_dict = json.load(open('./core/data/metanetx_sji.json','rb'))
    project = pd.read_excel(open(name+'_input.xlsx','rb'),sheet_name='project',header=None)
    project = project[0].values.tolist()[0]
    substrates = pd.read_excel(open(name+'_input.xlsx','rb'),sheet_name='substrates',index_col=0,header=None)
    substrates = substrates.to_dict()[1]
    products = pd.read_excel(open(name+'_input.xlsx','rb'),sheet_name='products',index_col=0,header=None)
    products = products.to_dict()[1]
    unknown_met = pd.read_excel(open(name+'_input.xlsx','rb'),sheet_name='novel_mets',header=None,names=['ID','ChEBI','SMILES'])
    if not unknown_met.empty:
        novel_mets = unknown_met['ID'].values.tolist()
    else:
        novel_mets = []
    primary_substrate = pd.read_excel(open(name+'_input.xlsx','rb'),sheet_name='primary_substrate',header=None)
    primary_substrate = primary_substrate[0].values.tolist()
    primary_product = pd.read_excel(open(name+'_input.xlsx','rb'),sheet_name='primary_product',header=None)
    primary_product = primary_product[0].values.tolist()
    product_name = pd.read_excel(open(name+'_input.xlsx','rb'),sheet_name='product_name',header=None)
    product_name = product_name[0].values.tolist()
    var_data = pd.read_excel(open(name+'_input.xlsx','rb'),sheet_name='var_data',header=None,index_col=0)
    iterations = var_data.loc['iterations'][1]
    distance = var_data.loc['distance'][1]
    nRxn1 = 0
    nRule = var_data.loc['nRule'][1]
    total_steps = nRule
    if len(novel_mets) != 0:
        get_molsig_exmetab(unknown_met)
    #pdb.set_trace()
    (solutions,nRxn2,v_maximum) = novoStoic_rule(reaction_dict,distance,substrates,products,novel_mets,
                iterations,project,nRxn1,nRule,total_steps)
    #pdb.set_trace()
    (met_result,sol_result) = intermediate_tracer(substrates,products,novel_mets,solutions,
                distance,reaction_dict, primary_substrate, primary_product,
                project)

    (solutions2,r_type,v_values) = novoStoic_rxn(reaction_dict,distance,substrates,products,novel_mets,
                iterations,project,nRxn2,nRxn2,nRxn2,met_result,primary_substrate,
                primary_product,solutions,v_maximum)

    (met_result,sol_result) = intermediate_tracer(substrates,products,novel_mets,solutions2,
                distance,reaction_dict, primary_substrate, primary_product,
                project)

    metab_dict = dict()
    
    for i in met_result:

        metab_dict[i] = met_iden(met_result[i],primary_substrate,primary_product)

    (r_type_input) = generate_files(sol_result,metab_dict,v_values,reaction_dict,r_type,product_name[0])
    out_ind = 0

    for i in sol_result:
        (out_ind) = structure_gen(sol_result[i],primary_substrate,solutions2[i],i,r_type[i],metab_dict[i],project,reaction_dict,r_type_input[i],product_name[0],out_ind)
        out_ind = out_ind + 1

def generate_files(sol_result,metab_dict,v_values,reaction_dict,r_type,product_name):
    output_index = 0
    num_rxn = pd.DataFrame(columns=['rxn','n_of_rxn'])
    r_type_input = dict()
    for i in sol_result:
        r_type_input[i] = dict()
        for j in sol_result[i]:
            r_type_input[i][j] = dict()
            reaction_output = pd.DataFrame(columns=['order', 'type'])
            rule = 0
            step = 0
            for k in range(0,len(sol_result[i][j])):
                if v_values[i][sol_result[i][j][k]] < 0:
                    flux_out = -1
                elif v_values[i][sol_result[i][j][k]] > 0:
                    flux_out = 1
                if metab_dict[i][j][k] in reaction_dict[sol_result[i][j][k]].keys() and r_type[i][sol_result[i][j][k]] == 'rxn':   
                    reaction_output2 = pd.DataFrame({'order':[k],'MNX_ID':[sol_result[i][j][k]],'type':['rxn'],'flux':[flux_out]})
                    r_type_input[i][j][k] = 'rxn'
                    step = step + 1

                else:
                    reaction_output2 = pd.DataFrame({'order':[k],'MNX_ID':[sol_result[i][j][k]],'type':['rule'],'flux':[flux_out]})
                    r_type_input[i][j][k] = 'rule'
                    rule = rule + 1
                    step = step + 1
                reaction_output = reaction_output.append(reaction_output2)
            reaction_output2 = pd.DataFrame({'order':[''],'MNX_ID':[''],'type':[''],'flux':['']})
            reaction_output = reaction_output.append(reaction_output2)
            reaction_output = reaction_output.sort_values(by=['order'])
            reaction_output2 = pd.DataFrame({'order':['metabolite'],'MNX_ID':['order'],'type':[''],'flux':['']})
            reaction_output = reaction_output.append(reaction_output2)
            for k in range(0,len(metab_dict[i][j])):
                reaction_output2 = pd.DataFrame({'order':[metab_dict[i][j][k]],'MNX_ID':[k],'type':[''],'flux':['']})
                reaction_output = reaction_output.append(reaction_output2)
            reaction_output.to_csv('./test/result/solution_'+product_name+'_'+str(step)+'_'+str(rule)+'_'+str(output_index)+'.csv',index=False)
            nr_input = ''.join('solutions_'+str(output_index))
            summation = sum(1 for x in r_type[i].values() if x == 'rxn')
            num_rxn2 = pd.DataFrame({'rxn': [nr_input],'n_of_rxn':[summation]})
            num_rxn = num_rxn.append(num_rxn2)
            output_index = output_index + 1
    num_rxn = num_rxn.sort_values('n_of_rxn',ascending=False)
    num_rxn.to_csv('./test/result/solution_order.csv',index=False)
    return r_type_input



    

def novoStoic_rule(model,radius,substrates,products,novel_mets,iterations,
                project,nRxn,nRule,nRxnRule):
    
    counter = 0
    solutions = dict()
    cofactors = pd.read_csv('./core/data/cofactors.csv',index_col=0)

    M = 4
    molsigs = pd.read_csv(os.path.join(data_dir,'molecular_signature_' + str(radius) + '.csv'),\
        index_col=0).fillna(0)
    rules = pd.read_csv(os.path.join(data_dir,'rule_' + str(radius) + '_noduplic.csv'),index_col=0)
    #reactions_index = list(pd.read_csv(os.path.join(data_dir,'rule_'+str(radius)+'.csv')))
    #metabolites_index = list(molsigs)
    ###### sets ############
    moiety_index = rules.index.tolist()     # moiety sets
    rl_list = list(pd.read_csv('./core//data/rule_0.csv', index_col=0).columns.values)
    for i in model.keys():
        if i not in rl_list: 
            del model[i]
    reactions_index = model.keys()
    #pd.read_csv(os.path.join(data_dir,'MNX_reactions.txt'),header=None)[0].tolist()
    rules_index = rules.columns.values.tolist()
    exchange_index = substrates.keys() + products.keys()

    ###### parameters ######
    # S(i,j) stoichiometry for each reaction
    S = get_S(model)
    metabolites_index = S.keys()
    metabolites_index.extend(novel_mets) # add novel metabolites not in DB
    #metabolites_index = list(pd.read_csv('./core/data/molecular_signature_0.csv', index_col=0).columns.values)
    #reactions_index = list(pd.read_csv('./core//data/rule_0_noduplic.csv', index_col=0).columns.values)
    #pdb.set_trace()
    #metabolites_index = [x.encode('UTF8') for x in metabolites_index]
    # T(m,r) contains atom stoichiometry for each rule
    T = rules.to_dict(orient='index')

    # C(m,i) contains moiety cardinality for each metaboli.to_dict(orient='index')
    C = molsigs.to_dict(orient='index')
    for m in moiety_index:
        C[m]['MNXM1'] = 0
        C[m]['MNXM195'] = 0

    for met in novel_mets:
        molsigs_product = pd.read_csv(project + '/novel_mets/molsig_'+ met + '_' + \
            str(radius) + '.csv', index_col=0)
        molsigs_product_dict = molsigs_product.to_dict(orient='index')
        for m in moiety_index:
            if m in molsigs_product_dict.keys():
                C[m][met] = molsigs_product_dict[m][met]
            else:
                C[m][met] = 0       

    ###### variables ######
    v = pulp.LpVariable.dicts("v", reactions_index, lowBound=-M, upBound=M,\
     cat='Integer')
    v_imb = pulp.LpVariable.dicts("v_imb", metabolites_index, lowBound=-M,\
     upBound=M, cat='Integer')
    v_rule = pulp.LpVariable.dicts("v_rule", rules_index, lowBound=-M,\
     upBound=M, cat='Integer')
    v_1 = pulp.LpVariable.dicts("v_1", rules_index, lowBound=0,\
     upBound=M, cat='Integer')
    v_2 = pulp.LpVariable.dicts("v_2", reactions_index, lowBound=0,\
     upBound=M, cat='Integer')
    v_EX = pulp.LpVariable.dicts("v_EX", exchange_index, lowBound=-M,\
     upBound=M, cat='Integer')
    # v_EX = pulp.LpVariable.dicts("v_EX", metabolites_index, lowBound=-M,\
    #  upBound=M, cat='Integer')
    y_rxn = pulp.LpVariable.dicts("y_rxn", reactions_index, lowBound=0,\
     upBound=1, cat='Binary')
    y_rule = pulp.LpVariable.dicts("y_rule", rules_index, lowBound=0,\
     upBound=1, cat='Binary')
    y_imb = pulp.LpVariable.dicts("y_imb", metabolites_index, lowBound=0,\
     upBound=1, cat='Binary')

    # create MILP problem
    lp_prob = pulp.LpProblem("novoStoic", pulp.LpMinimize)

    ####### objective ####
    #lp_prob += pulp.lpSum([v_1[r] for r in rules_index]) \
    #        + pulp.lpSum([v_2[r] for r in reactions_index]) \
            # + 10*pulp.lpSum([y_imb[i] for i in metabolites_index]) # test min imbalance
    lp_prob += pulp.lpSum([y_rxn[r] for r in reactions_index]) \
            + pulp.lpSum([y_rule[r] for r in rules_index]) \
    # lp_prob += pulp.lpSum([y_imb[i] for i in metabolites_index]), "MinRule"
    #new mass balance constraint

        ####### constraint ####
    # constraint 1: component balance
    for i in metabolites_index:
        if i in novel_mets:
            lp_prob += v_imb[i] + v_EX[i] == 0, 'mass_balance_novel_' + i
        elif i in exchange_index:
            lp_prob += pulp.lpSum([S[i][j] * v[j] for j in S[i].keys()]) ==\
                        v_imb[i] + v_EX[i],'mass_balance_' + i
        else:
            lp_prob += pulp.lpSum([S[i][j] * v[j] for j in S[i].keys()]) ==\
                        v_imb[i], 'mass_balance_' + i        

    # constraint 2: moiety change balance
    #pdb.set_trace()
    for m in moiety_index:
        lp_prob += pulp.lpSum([T[m][r] * v_rule[r] for r in rules_index]) + \
                   pulp.lpSum([C[m][i] * v_imb[i] for i in metabolites_index])==0,\
                   'moiety_balance_' + str(m)

    # constraint 3: constraint for exchange reactions
    for i,stoic in substrates.iteritems():
        if i != 'MNXM1':
            lp_prob += v_EX[i] == stoic, 'uptake' + i

    for i,stoic in products.iteritems():
        if i != 'MNXM1':
            lp_prob += v_EX[i] == stoic, 'produce' + i

    # for i in metabolites_index:
    #     if i not in exchange_index:
    #         lp_prob += v_EX[i] == 0, 'no_exchange' + i

    # constraint 4: number of reactions and rules
    LB = {}
    UB = {}
    for j in reactions_index:
        LB[j] = -M
        UB[j] = M

    for j in reactions_index:
            lp_prob += v[j] >= y_rxn[j]*LB[j]/2, "cons1_%s"%j
            lp_prob += v[j] <= y_rxn[j]*UB[j]/2, "cons2_%s"%j
    lp_prob += pulp.lpSum([y_rxn[j] for j in reactions_index]) <= nRxn, 'number_rxn'


    for j in rules_index:
            lp_prob += v_rule[j] >= y_rule[j]*LB[j], "cons1_rule_%s"%j
            lp_prob += v_rule[j] <= y_rule[j]*UB[j], "cons2_rule_%s"%j

    for j in rules_index:
            lp_prob += v_rule[j]  <= v_1[j], "abs_rule1_%s"%j
            lp_prob += -v_rule[j]  <= v_1[j], "abs_rule2_%s"%j
    
    for j in reactions_index:
            lp_prob += v[j]  <= v_2[j], "abs_reaction1_%s"%j
            lp_prob += -v[j] <= v_2[j], "abs_reaction2_%s"%j
            
    lp_prob += pulp.lpSum([y_rule[r] for r in rules_index]) <= nRule, 'number_rule'

    lp_prob += pulp.lpSum([y_rxn[j] for j in reactions_index]) + \
               pulp.lpSum([y_rule[r] for r in rules_index]) <= nRxnRule, 'total_number'

    # for i in metabolites_index:
    #         lp_prob += v_imb[i] >= y_imb[i]*(-M), "cons1_imb_%s"%i
    #         lp_prob += v_imb[i] <= y_imb[i]*M, "cons2_imb_%s"%i


    # R00479

    # lp_prob += pulp.lpSum([y_imb[r] for r in metabolites_index]) >= 5, 'number_imbmet'

    # lp_prob += pulp.lpSum([y_rxn[r] for r in reactions_index]) \
    #         + 10*pulp.lpSum([y_rule[r] for r in rules_index]) >=10, "at_least_one_rule"


    # constraint 5: thermodynamic feasibility

    ### solve 
    # pulp_solver = pulp.solvers.CPLEX_CMD(path=None, keepFiles=0, mip=1, msg=1,\
    #  options=['mip tolerances mipgap 0', 'mip tolerances absmipgap 0',\
    #   'mip tolerances integrality 0', 'simplex tolerances optimality 1E-9',\
    #   'simplex tolerances feasibility 1E-9',], timelimit=1200)
    # pulp_solver = pulp.solvers.GUROBI()

    pulp_solver = pulp.solvers.CPLEX_PY()
    #lp_prob.writeLP('kegg_novostoic_isobutanol.lp')
    
    for sol_num in range(0,iterations):
        integer_cut_rxns = []
        integer_cut_rules = []
        lp_prob.solve(pulp_solver)
        print("Status:", pulp.LpStatus[lp_prob.status], sol_num)

        if pulp.LpStatus[lp_prob.status] != 'Optimal': break

        print "----------sol: known rxns--------------",sol_num
        for j in reactions_index:
            if y_rxn[j].varValue >= 0.8:
                print j
                integer_cut_rxns.append(j)
        # pdb.set_trace()
        print "----------sol: rules--------------",sol_num
        for r in rules_index:
            if y_rule[r].varValue >= 0.8:
                integer_cut_rules.append(r)
                print r

        # T_new_all= postProcessRule(integer_cut_rules,model,cofactor,T,C,S,moiety_index)
        # get_pathway_from_sol(v_rule,v_EX,y_rule,integer_cut_rules,\
        #                 T_new_all,C,substrates_index,moiety_index)

        solutions[counter] = dict()
        count = 0
        nRxn2 = 0
        #with open(project + '/result/'+'/sol_fluxes_'+str(sol_num)+'.csv','wb') as f:
            #f.write('fluxes: \n')
        for j in reactions_index:
            if v[j].varValue != 0:
                #f.write(j + ','+ str(v[j].varValue) + ',' + str(y_rxn[j].varValue))
                #f.write('\n')
                solutions[counter][j] = v[j].varValue
                nRxn2 = nRxn2 + abs(v[j].varValue)
                count = count + 1
        #f.write('rules: \n')
        for r in rules_index:
            if v_rule[r].varValue !=0:
                #f.write(r + ','+ str(v_rule[r].varValue) + ',' + str(y_rule[r].varValue))
                #f.write('\n')
                solutions[counter][r] = v_rule[r].varValue
                nRxn2 = nRxn2 + abs(v_rule[r].varValue)
                count = count + 1
        #f.write('imbalance: \n')
        #for i in metabolites_index:
            #if v_imb[i].varValue !=0:
                #f.write(i + ','+ str(v_imb[i].varValue))
                #f.write('\n')
        #f.write('exchange:\n')
        #for j in exchange_index:
            #f.write(j + ','+ str(v_EX[j].varValue))
            #f.write('\n')
                
        length = len(integer_cut_rules) + len(integer_cut_rxns) - 1
        lp_prob += pulp.lpSum([y_rxn[r] for r in integer_cut_rxns ]) + \
                pulp.lpSum([y_rule[r] for r in integer_cut_rules ]) <= length,\
                             'intger_cult_' + str(counter)
        counter = counter + 1
    v_maximum = dict()
    for i in solutions.keys():
        v_maximum[i] = 0
        for j in solutions[i].keys():
            v_maximum[i] = v_maximum[i] + abs(solutions[i][j])


    return (solutions, nRxn2, v_maximum)

def novoStoic_rxn(model,radius,substrates,products,novel_mets,iterations,
                project,nRxn,nRule,nRxnRule,met_result,primary_substrate,
                primary_product,solutions,v_maximum):
    #pdb.set_trace()
    counter = 0
    solutions2 = dict()
    r_type = dict()
    flux = dict()
    sol_count = 1
    iterate = 0
    for things in met_result:
        #pdb.set_trace()
        true_metabolites = dict()
        
        for things2 in met_result[things]:
            for things3 in met_result[things][things2]:
                #for things4 in met_result[things][things2][things3]:
                true_metabolites[iterate] = met_result[things][things2][things3]
                iterate = iterate + 1
        #pdb.set_trace()
        sol_iter = solutions[things]
        #true_metabolites = json.load(open('intermediate_metabolites.json','rb'))
        cofactors = pd.read_csv('./core/data/cofactors.csv',index_col=0)
        M = 4
        molsigs = pd.read_csv(os.path.join(data_dir,'molecular_signature_' + str(radius) + '.csv'),\
            index_col=0).fillna(0)

        rules = pd.read_csv(os.path.join(data_dir,'rule_' + str(radius) + '_noduplic.csv'),index_col=0)
        #reactions_index = list(pd.read_csv(os.path.join(data_dir,'rule_'+str(radius)+'.csv')))
        #metabolites_index = list(molsigs)
        #metabolites_index.extend(novel_mets) # add novel metabolites not in DB
        #pdb.set_trace()
        ###### sets ############
        moiety_index = rules.index.tolist()     # moiety sets
        rl_list = list(pd.read_csv('./core//data/rule_0.csv', index_col=0).columns.values)
        for i in model.keys():
            if i not in rl_list: 
                del model[i]
        reactions_index = model.keys()
        rules_index = rules.columns.values.tolist()
        exchange_index = substrates.keys() + products.keys()

        #pdb.set_trace()
        for i in rules.columns.values.tolist():
            if i not in sol_iter:
                #pdb.set_trace()
                rules = rules.drop(i,axis = 1)
                rules_index.remove(i)


        #pdb.set_trace()

        ###### parameters ######
        # S(i,j) stoichiometry for each reaction
        S = get_S(model)
        metabolites_index = S.keys()
        metabolites_index.extend(novel_mets) # add novel metabolites not in DB
        

        # T(m,r) contains atom stoichiometry for each rule
        T = rules.to_dict(orient='index')

        # C(m,i) contains moiety cardinality for each metabolite
        C = molsigs.to_dict(orient='index')
        for m in moiety_index:
            C[m]['MNXM1'] = 0
            C[m]['MNXM195'] = 0

        for met in novel_mets:
            molsigs_product = pd.read_csv(project + '/novel_mets/molsig_'+ met + '_' + \
                str(radius) + '.csv', index_col=0)
            molsigs_product_dict = molsigs_product.to_dict(orient='index')
            for m in moiety_index:
                if m in molsigs_product_dict.keys():
                    C[m][met] = molsigs_product_dict[m][met]
                else:
                    C[m][met] = 0       

        ###### variables ######
        v = pulp.LpVariable.dicts("v", reactions_index, lowBound=-M, upBound=M,\
        cat='Integer')
        v_imb = pulp.LpVariable.dicts("v_imb", metabolites_index, lowBound=-M,\
        upBound=M, cat='Integer')
        v_rule = pulp.LpVariable.dicts("v_rule", rules_index, lowBound=-M,\
        upBound=M, cat='Integer')
        v_1 = pulp.LpVariable.dicts("v_1", rules_index, lowBound=0,\
        upBound=M, cat='Integer')
        v_2 = pulp.LpVariable.dicts("v_2", reactions_index, lowBound=0,\
        upBound=M, cat='Integer')
        v_EX = pulp.LpVariable.dicts("v_EX", exchange_index, lowBound=-M,\
        upBound=M, cat='Integer')
        # v_EX = pulp.LpVariable.dicts("v_EX", metabolites_index, lowBound=-M,\
        #  upBound=M, cat='Integer')
        y_rxn = pulp.LpVariable.dicts("y_rxn", reactions_index, lowBound=0,\
        upBound=1, cat='Binary')
        y_rule = pulp.LpVariable.dicts("y_rule", rules_index, lowBound=0,\
        upBound=1, cat='Binary')
        y_imb = pulp.LpVariable.dicts("y_imb", metabolites_index, lowBound=0,\
        upBound=1, cat='Binary')

        # create MILP problem
        lp_prob = pulp.LpProblem("novoStoic", pulp.LpMaximize)

        ####### objective ####
        lp_prob += pulp.lpSum([y_rxn[r] for r in reactions_index]) \
                # + 10*pulp.lpSum([y_imb[i] for i in metabolites_index]) # test min imbalance
        
        # lp_prob += pulp.lpSum([y_imb[i] for i in metabolites_index]), "MinRule"
        #new mass balance constraint
        met_flag = dict()
        for i in metabolites_index:
            met_flag[i] = 1

        metab_check = dict()
        
        #v_imb mass balance constraints
        for ind in metabolites_index:
            metab_check[ind] = dict()
            for ind2 in true_metabolites.keys():
                #pdb.set_trace()
                metab_check[ind][ind2] = 0
                for m in moiety_index:
                    metab_check[ind][ind2] = metab_check[ind][ind2] + abs(C[m][ind] - true_metabolites[ind2][m])
                if metab_check[ind][ind2] == 0:
                    met_flag[ind] = 0
        for ind in substrates:
            met_flag[ind] = 0
        for ind in products:
            met_flag[ind] = 0

        for i in metabolites_index:
            if met_flag[i] == 1:
                lp_prob += v_imb[i] == 0
        #pdb.set_trace()
            ####### constraint ####
        # constraint 1: component balance
        for i in metabolites_index:
            if i in novel_mets:
                lp_prob += v_imb[i] + v_EX[i] == 0, 'mass_balance_novel_' + i
            elif i in exchange_index:
                lp_prob += pulp.lpSum([S[i][j] * v[j] for j in S[i].keys()]) ==\
                            v_imb[i] + v_EX[i],'mass_balance_' + i
            else:
                lp_prob += pulp.lpSum([S[i][j] * v[j] for j in S[i].keys()]) ==\
                            v_imb[i], 'mass_balance_' + i        

        # constraint 2: moiety change balance

        for m in moiety_index:
            lp_prob += pulp.lpSum([T[m][r] * v_rule[r] for r in rules_index]) + \
                       pulp.lpSum([C[m][i] * v_imb[i] for i in metabolites_index])==0,\
                       'moiety_balance_' + str(m)


        # constraint 3: constraint for exchange reactions
        for i,stoic in substrates.iteritems():
            if i != 'MNXM1':
                lp_prob += v_EX[i] == stoic, 'uptake' + i

        for i,stoic in products.iteritems():
            if i != 'MNXM1':
                lp_prob += v_EX[i] == stoic, 'produce' + i

        # for i in metabolites_index:
        #     if i not in exchange_index:
        #         lp_prob += v_EX[i] == 0, 'no_exchange' + i

        # constraint 4: number of reactions and rules
        LB = {}
        UB = {}
        for j in reactions_index:
            LB[j] = -M
            UB[j] = M

        for j in reactions_index:
                lp_prob += v[j] >= y_rxn[j]*LB[j]/2, "cons1_%s"%j
                lp_prob += v[j] <= y_rxn[j]*UB[j]/2, "cons2_%s"%j
        lp_prob += pulp.lpSum([y_rxn[j] for j in reactions_index]) <= nRxn, 'number_rxn'
        lp_prob += pulp.lpSum([v_1[j] for j in rules_index]) + \
                   pulp.lpSum([v_2[r] for r in reactions_index]) <= v_maximum[things], 'maximum_steps'

        for j in rules_index:
                lp_prob += v_rule[j] >= y_rule[j]*LB[j], "cons1_rule_%s"%j
                lp_prob += v_rule[j] <= y_rule[j]*UB[j], "cons2_rule_%s"%j

        for j in rules_index:
            lp_prob += v_rule[j]  <= v_1[j], "abs_rule1_%s"%j
            lp_prob += -v_rule[j]  <= v_1[j], "abs_rule2_%s"%j
    
        for j in reactions_index:
            lp_prob += v[j]  <= v_2[j], "abs_reaction1_%s"%j
            lp_prob += -v[j] <= v_2[j], "abs_reaction2_%s"%j


        lp_prob += pulp.lpSum([y_rule[r] for r in rules_index]) <= nRule, 'number_rule'

        lp_prob += pulp.lpSum([y_rxn[j] for j in reactions_index]) + \
                   pulp.lpSum([y_rule[r] for r in rules_index]) <= nRxnRule, 'total_number'
        lp_prob += pulp.lpSum([y_rule[r] for r in rules_index]) >= 1, 'ensure valid solution'
        # for i in metabolites_index:
        #         lp_prob += v_imb[i] >= y_imb[i]*(-M), "cons1_imb_%s"%i
        #         lp_prob += v_imb[i] <= y_imb[i]*M, "cons2_imb_%s"%i

        # lp_prob +=v_rule['R01644'] ==1
        # lp_prob +=v_rule['R01486'] ==2

        # R00479

        # lp_prob += pulp.lpSum([y_imb[r] for r in metabolites_index]) >= 5, 'number_imbmet'

        # lp_prob += pulp.lpSum([y_rxn[r] for r in reactions_index]) \
        #         + 10*pulp.lpSum([y_rule[r] for r in rules_index]) >=10, "at_least_one_rule"


        # constraint 5: thermodynamic feasibility

        ### solve 
        # pulp_solver = pulp.solvers.CPLEX_CMD(path=None, keepFiles=0, mip=1, msg=1,\
        #  options=['mip tolerances mipgap 0', 'mip tolerances absmipgap 0',\
        #   'mip tolerances integrality 0', 'simplex tolerances optimality 1E-9',\
        #   'simplex tolerances feasibility 1E-9',], timelimit=1200)
        # pulp_solver = pulp.solvers.GUROBI()

        pulp_solver = pulp.solvers.CPLEX_PY()
        # lp_prob.writeLP('kegg_novostoic_isobutanol.lp')
    	y_store_rxns[things] = dict()
    	y_store_rules[things] = ditc()
        for sol_num in range(0,iterations):
            integer_cut_rxns = []
            integer_cut_rules = []
            lp_prob.solve(pulp_solver)
            print("Status:", pulp.LpStatus[lp_prob.status], sol_num)

            if pulp.LpStatus[lp_prob.status] != 'Optimal': break

            print "----------sol: known rxns--------------",sol_num
            for j in reactions_index:
                if v[j].varValue >= 0.8 or v[j].varValue <= -0.8:
                    print j
                    integer_cut_rxns.append(j)
            # pdb.set_trace()
            print "----------sol: rules--------------",sol_num
            for r in rules_index:
                if v_rule[r].varValue >= 0.8 or v_rule[r].varValue <= -0.8:
                    integer_cut_rules.append(r)
                    print r

            # T_new_all= postProcessRule(integer_cut_rules,model,cofactor,T,C,S,moiety_index)
            # get_pathway_from_sol(v_rule,v_EX,y_rule,integer_cut_rules,\
            #                 T_new_all,C,substrates_index,moiety_index)
            solutions2[counter] = dict()
            flux[counter] = dict()
            r_type[counter] = dict()
            #with open(project + '/result/'+'/sol_fluxes_'+str(sol_count)+'.csv','wb') as f:
            #    f.write('fluxes: \n')
            for j in reactions_index:
                if v[j].varValue != 0:
                    flux[counter][j] = v[j].varValue
            #            f.write(j + ','+ str(v[j].varValue) + ',' + str(y_rxn[j].varValue))
            #            f.write('\n')
                    solutions2[counter][j] = v[j].varValue
                    r_type[counter][j] = 'rxn'
            #    f.write('rules: \n')
            for r in rules_index:
                if v_rule[r].varValue !=0:
                    flux[counter][r] = v_rule[r].varValue
                        #f.write(r + ','+ str(v_rule[r].varValue) + ',' + str(y_rule[r].varValue))
                        #f.write('\n')
                    solutions2[counter][r] = v_rule[r].varValue
                    r_type[counter][r] = 'rule'
            #    f.write('imbalance: \n')
            #    for i in metabolites_index:
            #        if v_imb[i].varValue !=0:
            #            f.write(i + ','+ str(v_imb[i].varValue))
            #            f.write('\n')
            #    f.write('exchange:\n')
            #    for j in exchange_index:
            #        f.write(j + ','+ str(v_EX[j].varValue))
            #        f.write('\n')
            pdb.set_trace()
            y_store_rxn[things][sol_num] = integer_cut_rxns
            y_store_rules[thngs][sol_num] = integer_cut_rules    
            length = len(integer_cut_rules) + len(integer_cut_rxns) - 1
            lp_prob += pulp.lpSum([y_rxn[r] for r in integer_cut_rxns ]) + \
                    pulp.lpSum([y_rule[r] for r in integer_cut_rules ]) <= length,\
                                 'intger_cult_' + str(counter)
            
            counter = counter + 1
            sol_count = sol_count + 1
    return (solutions2,r_type,flux)
def postProcessRule(integer_cut_rules,model,cofactor,T,C,S,moiety_index):
    T_new_all = dict()
    T_new = dict()
    for r in integer_cut_rules:
        for m in moiety_index:
            T_new[m] = T[m][r]

        for i in cofactor:
            if S[i].get(r) is not None:
                for key in T_new.keys():
                    T_new[key] -= S[i][r]*C[key][i]
        T_new_all[r] = T_new
    return T_new_all

def get_pathway_from_sol(v_rule,v_EX,y_rule,integer_cut_rules,\
                        T_new_all,C,substrates_index,moiety_index):
    startMoieties = dict()
    for m in moiety_index:
        startMoieties[m] = 0
    for i in substrates_index:
        for m in moiety_index:
            startMoieties[m] += C[m][i]*v_EX[i]

    length = len(integer_cut_rules)
    # pdb.set_trace()
    for pathway in itertools.permutations(integer_cut_rules,length):
        indicator = 1
        for ridx in pathway:
            for m in moiety_index:
                startMoieties[m] += T_new_all[ridx][m]*v_rule[ridx]
                if startMoieties[m] < 0:
                    indicator = 0
                    break
        if indicator == 1:
            print 'feasible pathway', pathway



def map_to_KEGG(model):
    metList = [ met for met in model.metabolites]
    # print metList
    # pdb.set_trace()

    # convert to KEGG id
    for m in model.metabolites:
        # print m.id,':'
        if m.id == 'nh4_c' or m.id == 'nh4_e' or m.id == 'q8h2_c': continue
        
        annot = m.annotation['kegg.compound']
        # print m.id, annot
        if not isinstance(annot, basestring):
            m.id = annot[0]
        else:
            m.id = annot
        # print m.id
    metList = [ met.id for met in model.metabolites]
    # print metList
    model.repair()

    # print len(model.reactions)

    delete_exchange = ['BIOMASS_Ecoli_core_w_GAM']
    for r in model.reactions:
        if 'EX_' in r.id: delete_exchange.append(r.id)
    model.remove_reactions(delete_exchange)

    # print len(model.reactions)
    return model

def model_remove_ex(model):
    delete_exchange = ['BIOMASS_Ecoli_core_w_GAM']
    for r in model.reactions:
        if 'EX_' in r.id: delete_exchange.append(r.id)
    model.remove_reactions(delete_exchange)

    # print len(model.reactions)
    return model


if __name__ == '__main__':
    model = cobra.io.read_sbml_model('./e_coli_core.xml')

    model.optimize()
    print model.solution.status
    print model.solution.f

    model = model_remove_ex(model)

    novoStoic(model,2)
