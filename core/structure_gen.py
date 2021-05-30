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
data_dir = './core/data'
radius = 1
def structure_gen(rules,substrate,solution,iteration,r_type,metab_dict,project,reaction_dict,r_type_input,product_name,out_ind):
	filename = 'filename'
	#smarts = pd.read_csv(os.path.join(data_dir,'retrorule.csv'))
	#smiles = pd.read_csv(os.path.join(data_dir,'smiles.csv'))
	#met = smiles.iloc[smiles.iloc[:,0].tolist().index(substrate[0]),1]
	#reaction_dict = json.load(open('./core/data/optstoic_v3_reduced_Sji.json','rb'))
	#met = 'CC(C)(COP(=O)(O)OP(=O)(O)OCC1C(C(C(O1)N2C=NC3=C2N=CN=C3N)O)OP(=O)(O)O)C(C(=O)NCCC(=O)NCCSC(=O)CCC(=O)O)O'
	#rule_test = '([C:1]-[S:2]-[C:3](-[C:4])=[O:5])>>([C:1]-[S:2]-[H].[C:4]-[C:3](=[O:5])-[H])'
	#pdb.set_trace()
	r_type_check = dict()
	for j in rules.keys():
		#pdb.set_trace()
		r_type_check[j] = 0
		checker = 0
		store_rxn = dict()
		for k in range(0,len(rules[j])-1):
			store_rxn[rules[j][k]] = 0
		for k in range(0,len(rules[j])-1):
			if r_type[rules[j][k]] == 'rxn':
				store_rxn[rules[j][k]] = store_rxn[rules[j][k]] + 1
				if store_rxn[rules[j][k]] == 1:
					r_type_check[j] = r_type_check[j] + 1
				if metab_dict[j][k] in reaction_dict[rules[j][k]].keys():
					checker = checker + 1;
			#pdb.set_trace()
		if r_type_check[j] == checker:

			rule = dict()
			indices = dict()	
			for i in range(0,len(rules[j])):
				sep = '_'
				#pdb.set_trace()
				bb = rules[j][i]
				bb2 = bb.split(sep,1)[0]
				bb3 = solution[bb2]
				if bb3 > 0:
					bb3 = 1
				elif bb3 < 0:
					bb3 = -1
				#pdb.set_trace()
				#indices = [k for k, s in enumerate(smarts.iloc[:,0].tolist()) if bb2 in s]
				#pdb.set_trace()
				#for jj in indices:
				#	if smarts.iloc[jj,1] == bb3:
				#		rule[i] = AllChem.ReactionFromSmarts(smarts.iloc[jj,2])
			#if len(rule) == len(rules[j]):
			#	for i in range(0,len(rule)):
			#		if i == 0:
			#			#pdb.set_trace()
			#			met = Chem.MolFromSmiles(met)
			#			met = AddHs(met)
			#			Draw.MolToFile(met,'./metabolite_'+str(iteration)+'_'+str(j)+'_'+str(i)+'.png',size=(1500, 1500))
			#			images2 = Image.open('./metabolite_'+str(iteration)+'_'+str(j)+'_'+str(i)+'.png')
			#			bg = PIL.Image.new(images2.mode, images2.size, images2.getpixel((0, 0)))
			#			diff = PIL.ImageChops.difference(images2, bg)
			#			bbox = diff.getbbox()
			#			images2 = images2.crop(bbox).save('./metabolite_'+str(iteration)+'_'+str(j)+'_'+str(i)+'.png')
			#			product = rule[i].RunReactants((met,))[0][0]
			#			product = Chem.MolToSmiles(product)
			#			product = Chem.MolFromSmiles(product)
			#			product = AddHs(product)
			#			Draw.MolToFile(product,'./metabolite_'+str(iteration)+'_'+str(j)+'_'+str(i+1)+'.png',size=(1500, 1500))
			#			images2 = Image.open('./metabolite_'+str(iteration)+'_'+str(j)+'_'+str(i+1)+'.png')
			#			bg = PIL.Image.new(images2.mode, images2.size, images2.getpixel((0, 0)))
			#			diff = PIL.ImageChops.difference(images2, bg)
			#			bbox = diff.getbbox()
			#			images2 = images2.crop(bbox).save('./metabolite_'+str(iteration)+'_'+str(j)+'_'+str(i+1)+'.png')
			#			met = Chem.MolToSmiles(met)
			#			#pdb.set_trace()
			#		else:
			#			product = rule[i].RunReactants((product,))[0][0]
			#			product = Chem.MolToSmiles(product)
			#			product = Chem.MolFromSmiles(product)
			#			product = AddHs(product)
			#			Draw.MolToFile(product,'./metabolite_'+str(iteration)+'_'+str(j)+'_'+str(i+1)+'.png',size=(1500, 1500))
			#			images2 = Image.open('./metabolite_'+str(iteration)+'_'+str(j)+'_'+str(i+1)+'.png')
			#			bg = PIL.Image.new(images2.mode, images2.size, images2.getpixel((0, 0)))
			#			diff = PIL.ImageChops.difference(images2, bg)
			#			bbox = diff.getbbox()
			#			images2 = images2.crop(bbox).save('./metabolite_'+str(iteration)+'_'+str(j)+'_'+str(i+1)+'.png')
			#			#pdb.set_trace()
			#	#pdb.set_trace()
			#	full_image = 1
			#	image_construct(rules[j],metab_dict[j],solution,project,filename,iteration,project,j,full_image,r_type_input[j],out_ind,product_name)
			#	out_ind = out+ind + 1
			#else:
			full_image = 0
			image_construct(rules[j],metab_dict[j],solution,project,filename,iteration,project,j,full_image,r_type_input[j],out_ind,product_name)
			out_ind = out_ind + 1
	return out_ind

def image_construct(sol_result,met_result,solution,title,filename,iteration,project,i,full_image,r_type_input,out_ind,product_name):
	#met_names = met_iden(met_results)
	sep = '_'
	images = dict()
	#pdb.set_trace()
	color_configs = {}
	colorConfig = dict(COFACTOR_SHAPE="ellipse",  # "diamond"
          		                OTHER_COFACTOR_COLOR="#7F7F7F",
               		            NONCOFACTOR_SHAPE="plaintext",  # "box"
                      		    NONCOFACTOR_COLOR="transparent",  # "#CCFF33"
                           		REACTION_COLOR="black",
                           		EDGE_COLOR="black",  # "#505050"
                           		RXN_NODE_COLOR="black",
                           		BACKGROUND_COLOR="transparent",
                           		ALL_FONT_COLOR="white")
	if full_image == 1:
		for j in met_result:
			images2 = map(Image.open,['metabolite_'+str(iteration)+'_'+str(i)+'_'+str(j)+'.png'])
			imageFormat = 'png'
			g = gv.Digraph('G', format='png', engine='dot')
			g.graph_attr['rankdir'] = "LR"
			g.graph_attr['size'] = "10, 10"
			if imageFormat == 'png':
				g.graph_attr['dpi'] = '300'
			elif imageFormat == 'svg':
				g.graph_attr['dpi'] = '72'
			g.graph_attr['forcelabels'] = 'true'
			g.graph_attr['labelloc'] = 't'  # top or 'b'
			h = gv.Digraph('G', format='png', engine='dot')
			h.graph_attr['rankdir'] = "LR"
			h.graph_attr['size'] = "10, 10"
			if imageFormat == 'png':
				h.graph_attr['dpi'] = '300'
			elif imageFormat == 'svg':
				h.graph_attr['dpi'] = '72'
			h.graph_attr['forcelabels'] = 'true'
			h.graph_attr['labelloc'] = 't'  # top or 'b'

			d = gv.Digraph('G', format='png', engine='dot')
			d.graph_attr['rankdir'] = "BT"
			d.graph_attr['size'] = "10, 10"
			if imageFormat == 'png':
				d.graph_attr['dpi'] = '300'
			elif imageFormat == 'svg':
				d.graph_attr['dpi'] = '72'
			d.graph_attr['forcelabels'] = 'true'
			d.graph_attr['labelloc'] = 't'  # top or 'b'
			d.graph_attr['label'] = '1,4 butanediol synthesis'
			lineW = '3'

		for k in met_result.keys():
			g.node('metabolite_'+str(iteration)+'_'+str(i)+'_'+str(k),shape="none",label="",image='metabolite_'+str(iteration)+'_'+str(i)+'_'+str(k)+'.png')
		#pdb.set_trace()
		rule = 0
		step = 0
		for k in range(0,len(sol_result)):
			g.edge('metabolite_'+str(iteration)+'_'+str(i)+'_'+str(k),'metabolite_'+str(iteration)+'_'+str(i)+'_'+str(k+1), weight='1',
   	        	           	penwidth=lineW, color=colorConfig['EDGE_COLOR'],label=''.join(sol_result[k]+'\n'+r_type_input[k]))
			if rule_input_type[k] == 'rxn':
				step = step + 1
			else:
				step = step + 1
				rule = rule + 1
		g.render(filename, cleanup=True)
		#pdb.set_trace()
		for k in met_result.keys():
			h.node('metabolite_'+str(iteration)+'_'+str(i)+'_'+str(k),shape="ellipse",label=met_result[k])
		for k in range(0,len(sol_result)):
			h.edge('metabolite_'+str(iteration)+'_'+str(i)+'_'+str(k),'metabolite_'+str(iteration)+'_'+str(i)+'_'+str(int(k)+1), weight='1',
   	        	           	penwidth=lineW, color=colorConfig['EDGE_COLOR'],label=''.join(sol_result[k]+'\n'+r_type_input[k]))
		h.render('other',cleanup=True)
		d.node('image2',shape='none',label="",image='other.png')
		d.node('image1',shape="none",label="",image=filename+'.png')
		d.edge('image1','image2', weight='0',
   	    	               penwidth='0', arrowhead="none",color=colorConfig['EDGE_COLOR'],label="")
		d.render(project+'/image_'+product_name+'_'+str(step)+'_'+str(rule)+'_'+str(out_ind),cleanup=True)
	else:
		imageFormat = 'png'
		h = gv.Digraph('G', format='png', engine='dot')
		h.graph_attr['rankdir'] = "LR"
		h.graph_attr['size'] = "10, 10"
		if imageFormat == 'png':
			h.graph_attr['dpi'] = '300'
		elif imageFormat == 'svg':
			h.graph_attr['dpi'] = '72'
		h.graph_attr['forcelabels'] = 'true'
		h.graph_attr['labelloc'] = 't'  # top or 'b'
		h.graph_attr['label'] = '1,4 butanediol synthesis'
		lineW = '3'
		for k in met_result.keys():
			h.node('metabolite_'+str(iteration)+'_'+str(i)+'_'+str(k),shape="ellipse",label=met_result[k])
		step = 0
		rule = 0
		for k in range(0,len(sol_result)):
			h.edge('metabolite_'+str(iteration)+'_'+str(i)+'_'+str(k),'metabolite_'+str(iteration)+'_'+str(i)+'_'+str(int(k)+1), weight='1',
   	       		           	penwidth=lineW, color=colorConfig['EDGE_COLOR'],label=''.join(sol_result[k]+'\n'+r_type_input[k]))
			if r_type_input[k] == 'rxn':
				step = step + 1
			else:
				rule = rule + 1
				step = step + 1
		h.render(project+'/image_'+product_name+'_'+str(step)+'_'+str(rule)+'_'+str(out_ind),cleanup=True)

#for i in met_result.keys():
#	for j in met_result[i].keys():
#		os.remove('metabolite_'+str(iteration)+'_'+str(i)+'_'+str(j)+'.png')
	#os.remove(filename+'.png')
	#os.remove('other.png')
def met_iden(met_results,primary_substrate,primary_product):
	molsigs = pd.read_csv(os.path.join(data_dir,'molecular_signature_' + str(radius) + '.csv'),\
            index_col=0).fillna(0)
	moiety_index = molsigs.index.tolist()
	met_names = dict()
	gg = 0
	#pdb.set_trace()
	for k in met_results.keys():
		met_names[k] = dict()
		iterate = 0
		#pdb.set_trace()
		for i in met_results[k].keys():
			#pdb.set_trace()
			if i == 0:
				met_names[k][i] = primary_substrate[0]
			elif i == len(met_results[k].keys()) - 1:
				met_names[k][i] = primary_product[0]
			else:
				solfind = 0
				for j in molsigs:
					check_iden = 0
					for m in moiety_index:
						check_iden = check_iden + abs(met_results[k][i][m]-molsigs[j][m])
					if check_iden == 0:
						#pdb.set_trace()
						met_names[k][i] = j
						#iterate = iterate + 1
						solfind = 1
				if solfind == 0:
					met_names[k][i] = 'Hypothetical_Metabolite_'+str(iterate)
					iterate = iterate + 1
			gg = gg + 1
			#pdb.set_trace()
	#pdb.set_trace()
	return met_names




