
#Derived from B&B functions and some other home-made functions for QoL
#These scripts are for manipulating models easily

import os
import pandas as pd
import cobra
from datetime import datetime
import xlsxwriter




#Can be used to set fixed bounds or two bounds
def set_flux_bounds(r_id, bounds, model):
	if len(bounds) == 2:
		r_obj = model.reactions.get_by_id(r_id)
		r_obj.bounds = bounds
	else:
		print("Please include lower/upper bounds in tuple format")


#Generate constraint
def set_fix_flux_ratio(r_dict, model):
	if len(r_dict) == 2: #Needs to be 2
		key_list = list(r_dict.keys())
		value_list = list(r_dict.values())

		#First Reaction
		r0_id = key_list[0]
		r0_obj = model.reactions.get_by_id(r0_id)
		r0_flux = value_list[0]

		#Second Reaction
		r1_id = key_list[1]
		r1_obj = model.reactions.get_by_id(r1_id)
		r1_flux = value_list[1]

	#Adding ratio constraints to model
		constraint = model.problem.Constraint(r0_flux * r1_obj.flux_expression - r1_flux * r0_obj.flux_expression, lb=0, ub=0)
		#model.add_cons_vars(constraint)

		#Return constraint
		return constraint

#Medium needs to be formatted 
#Model is read from csv and needs to be formatted in the following format:
#exchange reaction | upper bound


def read_medium_csv(path, model):
	
	if os.path.exists(path):
		read = pd.read_csv(path, encoding='latin-1')
		medium_dict = dict(read.values)
		model.medium = medium_dict
		return model.medium
		
	else:
		print("Path doesn't exist")

	
def get_rxn(model, reaction_name):
	if model:
		if model.reactions.get_by_id(reaction_name):
			return model.reactions.get_by_id(reaction_name)
		else:
			print("No reaction exists")

	else:
		print("Model doesn't exist!")




def get_flux_exp(model, reaction_name):

	if model:
		if model.reactions.get_by_id(reaction_name.id):
			reaction_flux = model.reactions.get_by_id(reaction_name.id).flux_expression

			return reaction_flux
		else:
			print("No reaction exists")

	else:
		print("Model doesn't exist!")

def set_flux_to_1e6(model):
	if  model:
		for rxns in model.reactions:
		    if rxns.bounds != (0,0): #If reactions not 0,0
		        if rxns.upper_bound != 1e6 and rxns.upper_bound !=  0:
		            rxns.upper_bound =  1e6
		        if rxns.lower_bound != -1e6 and rxns.lower_bound != 0:
		            rxns.lower_bound =  -1e6


def get_mets(model, met_id):
	if model:
		if model.metabolites.get_by_id(met_id):
			return model.reactions.get_by_id(met_id)
		else:
			print("met doesn't exist")
	else:
		print("Model doesn't exist!")
		
def save_fba_matrix(filename, flux_matrix, path):		
	now = datetime.now()
	date_string = now.strftime("%Y%m%d-%H:%M")


	filename = '{}-{}'.format(filename, date_string)

	directory = os.path.dirname(path)

	if not os.path.exists(directory):
		os.makedirs(directory)
		print(str(directory))

	filepath = str(str(path) + str(filename) + '.tsv')
	file = open(filepath, 'w')
	flux_matrix.to_csv(file, sep='\t')
	file.close()
	print("Successfully saved {} to {}".format(filename, filepath))


def create_xlsx(filename, path):		
	now = datetime.now()
	date_string = now.strftime("%Y%m%d-%H:%M")


	filename = '{}-{}'.format(filename, date_string)

	directory = os.path.dirname(path)

	if not os.path.exists(directory):
		os.makedirs(directory)
		print(str(directory))

	filepath = str(str(path) + str(filename) + '.xlsx')
	workbook = pd.ExcelWriter(filepath, engine='xlsxwriter')

	print("Successfully generated {} in {}".format(filename, filepath))
	return workbook
			