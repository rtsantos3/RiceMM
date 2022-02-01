
#Derived from B&B functions and some other home-made functions for QoL
#These scripts are for manipulating models easily

import os
import pandas as pd


#Can be used to set fixed bounds or two bounds
def set_flux_bounds(r_id, bounds, model):
	if len(bounds) == 2:
		r_obj = model.reactions.get_by_id(r_id)
		r_obj.bounds = bounds
	else:
		print("Please include lower/upper bounds in tuple format")



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
		constraint = model.problem.Constraint(r0_flux * r1_obj.flux_expression - r1_flux * r0_obj.flux_expression, lb=0, ub=1000)
		model.add_cons_vars(constraint)

		#Return constraint
		return constraint

#Medium needs to be formatted 
#Model is read from csv and needs to be formatted in the following format:
#exchange reaction | upper bound


def read_medium_csv(path, model):
	
	if os.path.exists(path):
		read = pd.read_csv(path)
	medium_dict = dict(read.values)
	model.medium = medium_dict
	return model.medium

