
#Derived from B&B functions
#These scripts are for manipulating models easily


#Can be used to set fixed bounds or two bounds
def set_flux_bounds(r_id, bounds, model):
	if len(bounds) == 2:
		r_obj = model.reactions.get_by_id(r_id)
		r_obj.bounds = bounds
	else:
		print("Please include lower/upper bounds in tuple format")



def set_fix_flux_ratio(r_dict, model):
	if len(r_dict) == 2: #Needs to be 2
		#First Reaction
		r1_id = r_dict.keys()[0] 
		r1_obj = model.reactions.get_by_id(r1_id)
		r1_flux = r_dict.values()[0] 

		#Second Reaction
		r2_id = r_dic.keys()[1]
		r2_obj = model.reactions.get_by_id(r2_id)
		r2_flux = r_dict.values()[1]

		#Adding ratio constraints to model
		constraint = model.problem.constraint(r1_flux * r2_obj.flux_expression - r2_flux * r1_obj.flux_expression, lb=0, ub=0)
		model.add_cons_vars(constraint)

		#Return constraint
		return constraint

