#FUNCTIONS

###################################################
# This function estimates maintenance ATPase flux #
# based on Fraser's finding that maintenance      #
# increases linearly with PPFD                    #
###################################################
def estimateMaintenance(PPFD):
    ATPase = (PPFD*0.0049)+2.7851 #Equation form : y=mx+c
    return ATPase


####################################################
# This function updates Biomass composition of a   #
# model from a python dictionary                   #
####################################################
def updateBiomassComposition(model,BiomComp):
    newBiomass = dict()
    fin = open(BiomComp,"r")
    for line in fin:
        line=line.replace("\n","")
        lineparts = line.split(",")
        #print lineparts
        newBiomass[lineparts[0]]=float(lineparts[1])
    for i in [1,2]:
        temp = dict(model.reactions.get_by_id("Biomass_tx"+str(i)).metabolites)
        for met in temp.keys():
            model.reactions.get_by_id("Biomass_tx"+str(i)).add_metabolites({met:-1*temp.get(met)})
            #print(model.reactions.get_by_id("Biomass_tx"+str(i)).reaction)
        for met in newBiomass.keys():
            if met.startswith("S_"):
                continue
            if met.__contains__("TP") or met.__contains__("pTRP") or met.__contains__("pCYS") or met.__contains__("pPRO"):
                continue
                print met
            model.reactions.get_by_id("Biomass_tx"+str(i)).add_metabolites({model.metabolites.get_by_id(met+str(i)):-1*newBiomass.get(met)})
    return model


###################################################################
# This function writes a set of solutions and FVA results to file #
###################################################################
def printFluxesAndFVA(model,solutiondict,FVAdict,outfile):
    fout = open(outfile,"w")
    fout.write("Rxn ID\tReaction name\tEquation")
    for PPFD in sorted(solutiondict.keys()):
        fout.write("\tFlux, PPFD="+str(PPFD)+"\tMax Flux, PPFD="+str(PPFD)+"\tMin Flux, PPFD="+str(PPFD))
    fout.write("\n")
    for rxn in model.reactions:
        fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction)
        for PPFD in sorted(solutiondict.keys()):
            fout.write("\t"+str(solutiondict[PPFD][rxn.id])+"\t"+str(FVAdict[PPFD]["maximum"][rxn.id])+"\t"+str(FVAdict[PPFD]["minimum"][rxn.id]))
        fout.write("\n")
    fout.close()


######################################################
# This function estimates biomass/phloem output flux #
# at which the net CO2 uptake rate is equal to the   #
# user defined value                                 #
####################################################
def estimateOutputFromNetCO2(model,netCO2uptake,Output_ID="diel_biomass",Vc_ID="RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1",CO2in_ID="CO2_tx1",verbose=False):

    from cobra import flux_analysis
    # Initally constraint Vc flux to net CO2 uptake rate
    model.reactions.get_by_id(Vc_ID).lower_bound = netCO2uptake
    model.reactions.get_by_id(Vc_ID).upper_bound = netCO2uptake

    if verbose:
        print("Photons ="+str(model.reactions.get_by_id("Photon_tx1").upper_bound))
        print("ATPase ="+str(model.reactions.get_by_id("ATPase_tx1").upper_bound))
        print("Vc flux"+str(model.reactions.get_by_id(Vc_ID).upper_bound))

    #perform pFBA
    flux_analysis.parsimonious.pfba(model)

    if verbose:
        print("CO2 uptake ="+str(model.reactions.get_by_id("CO2_tx1").x))


    #unconstrain Vc
    model.reactions.get_by_id(Vc_ID).lower_bound = 0
    model.reactions.get_by_id(Vc_ID).upper_bound = 1000

    #set loop counter
    i=0

    #Use a while loop to increase Vc flux until net CO2 rate is similar to given value (or loop counter hits 10)
    while((netCO2uptake - model.reactions.get_by_id(CO2in_ID).x)/netCO2uptake > 0.0001 and i<10):
        i=i+1
        prev = model.reactions.get_by_id(Output_ID).x
        # Increment in Vc flux is set by given netCo2 uptake - model predicted CO2 uptake rate in previous pFBA run
        now = prev + (prev*((netCO2uptake - model.reactions.get_by_id(CO2in_ID).x)/netCO2uptake))
        model.reactions.get_by_id(Output_ID).lower_bound = now
        model.reactions.get_by_id(Output_ID).upper_bound = now

        flux_analysis.parsimonious.pfba(model)
        if verbose:
            print("----"+str(i)+"----")
            print("Vc flux ="+str(model.reactions.get_by_id(Vc_ID).x))
            print("net CO2 uptake ="+str(model.reactions.get_by_id(CO2in_ID).x))
            print("Target CO2 uptake ="+str(netCO2uptake))
            print("Before:"+str(prev))
            print("After:"+str(now))
        if i>10:
            print("Warning: Loop counter greater than 10")
    return prev


########################################################
#This function was used to set up a C3 leaf diel model #
########################################################
def setupDielModel(core_model,transferMets):
    from cobra.core import Metabolite, Reaction
    import re

    #create two copies of model elements for day and night
    cobra_model2 = core_model.copy()
    for met in cobra_model2.metabolites:
        met.id = met.id+"1"
        met.compartment = met.compartment+"1"
    for rxn in cobra_model2.reactions:
        rxn.id = rxn.id+"1"

    cobra_model3 = core_model.copy()
    for met in cobra_model3.metabolites:
        met.id = met.id+"2"
        met.compartment = met.compartment+"2"
    for rxn in cobra_model3.reactions:
        rxn.id = rxn.id+"2"

    #merge the day and night model
    cobra_model = cobra_model2+cobra_model3
    for met in cobra_model3.metabolites:
        if not cobra_model.metabolites.__contains__(met.id):
            cobra_model.add_metabolites(met.copy())

    met1 = Metabolite("X_Phloem_contribution_t1",name="Phloem output during the day")
    cobra_model.reactions.get_by_id("Phloem_output_tx1").add_metabolites({met1:1})
    met2 = Metabolite("X_Phloem_contribution_t2",name="Phloem output during at night")
    cobra_model.reactions.get_by_id("Phloem_output_tx2").add_metabolites({met2:1})

    rxn = Reaction("diel_biomass")
    rxn.add_metabolites({met1:-3,met2:-1})
    rxn.lower_bound = 0
    rxn.upper_bound = 1000
    cobra_model.add_reaction(rxn)

    #Adding reactions to allow for day-night metabolite accumulations
    tmfile = open(transferMets,"r")
    tmset=set()
    for line in tmfile:
        tmset.add(line.replace("\n",""))

    for met in tmset:
        if met == "AMMONIUM_v" or met=="FRUCTAN_v":
            continue
        tempRxn = Reaction(met+"_dielTransfer")
        tempRxn.add_metabolites({cobra_model.metabolites.get_by_id(met+"1"):-1,cobra_model.metabolites.get_by_id(met+"2"):1})
        tempRxn.lower_bound=-1000
        if not ((met == "STARCH_p") or (met == "SUCROSE_v") or (met == "MAL_v") or (met == "aMAL_v") or (met == "NITRATE_v") or (met == "CIT_v") or (met == "aCIT_v") or (met == "PROTON_v")):
            tempRxn.lower_bound=0
        tempRxn.upper_bound=1000
        cobra_model.add_reaction(tempRxn)

    fractionMets=dict()
    for rxn in cobra_model.reactions:
        for met in rxn.metabolites.keys():
            prefix=""
            a=re.search("^a{1,3}",met.id)
            anion=""
            if a:
                anion=a.group(0)
                prefix=anion
            b=re.search("^b{1,3}",met.id)
            basic=""
            if b:
                basic=b.group(0)
                prefix=basic
            if ((not prefix == "") and met.compartment == "v1"):
                fractionMets[met]=prefix

    temp=cobra_model.copy()
    for met in fractionMets.keys():
        for rxn in met.reactions:
            if rxn.id.__contains__("_dielTransfer"):
                continue
            else:
                mainMet = met.id[len(fractionMets[met]):]
                coeff1 = temp.reactions.get_by_id(rxn.id).metabolites.get(temp.metabolites.get_by_id(mainMet))
                coeff2 = temp.reactions.get_by_id(rxn.id).metabolites.get(temp.metabolites.get_by_id(met.id))
                if not coeff1:
                    coeff1=0
                if not coeff2:
                    coeff2=0
                total = coeff1 + coeff2
                coeff1 = float(coeff1)/total
                coeff2 = float(coeff2)/total
                if cobra_model.reactions.has_id(met.id[0:len(met.id)-1]+"_dielTransfer"):
                    ub = temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").upper_bound
                    lb = temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").lower_bound
                    temp.reactions.get_by_id(met.id[0:len(met.id)-1]+"_dielTransfer").remove_from_model()
                    temp.reactions.get_by_id(mainMet[0:len(mainMet)-1]+"_dielTransfer").remove_from_model()
                    Reac = Reaction(mainMet[0:len(mainMet)-1]+"_dielTransfer",name=mainMet+"_dielTransfer")
                    Reac.add_metabolites({temp.metabolites.get_by_id(met.id[0:len(met.id)-1]+"1"):-coeff2,temp.metabolites.get_by_id(met.id[0:len(met.id)-1]+"2"):coeff2,temp.metabolites.get_by_id(mainMet[0:len(mainMet)-1]+"1"):-coeff1,temp.metabolites.get_by_id(mainMet[0:len(mainMet)-1]+"2"):coeff1})
                    Reac.lower_bound=lb
                    Reac.upper_bound=ub
                    temp.add_reaction(Reac)
                    print Reac.reaction
                break
    ####ADD CONSTRAINTS TO MODEL####
    cobra_model = temp.copy()

    #objective function
    cobra_model.reactions.get_by_id("diel_biomass").objective_coefficient=1
    #Leaves - light
    cobra_model.reactions.get_by_id("Sucrose_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("Sucrose_tx1").upper_bound=0
    cobra_model.reactions.get_by_id("GLC_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("GLC_tx1").upper_bound=0
    cobra_model.reactions.get_by_id("CO2_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("NH4_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("NH4_tx1").upper_bound=0
    #Leaves - dark
    cobra_model.reactions.get_by_id("Sucrose_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("Sucrose_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("GLC_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("GLC_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("Photon_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("Photon_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("NH4_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("NH4_tx2").upper_bound=0
    cobra_model.reactions.get_by_id("CO2_tx2").upper_bound=0

    #Set pG6P transporter to 0
    cobra_model.reactions.get_by_id("G6P_Pi_pc1").lower_bound=0
    cobra_model.reactions.get_by_id("G6P_Pi_pc1").upper_bound=0
    cobra_model.reactions.get_by_id("G6P_Pi_pc2").lower_bound=0
    cobra_model.reactions.get_by_id("G6P_Pi_pc2").upper_bound=0

    #Turn off PTOX
    cobra_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").lower_bound=0
    cobra_model.reactions.get_by_id("Plastoquinol_Oxidase_p1").upper_bound=0

    #nitrate uptake constrain
    Nitrate_balance = Metabolite("Nitrate_bal_c", name = "Weights to balance nitrate uptake", compartment = "c1")
    cobra_model.reactions.get_by_id("Nitrate_ec1").add_metabolites({Nitrate_balance:-2})
    cobra_model.reactions.get_by_id("Nitrate_ec2").add_metabolites({Nitrate_balance:3})

    #Rubisco balance
    Rubisco_balance = Metabolite("rubisco_bal_p1", name = "Weights to balance RuBP carboxygenase oxygenase balance", compartment = "p1")
    cobra_model.reactions.get_by_id("RXN_961_p1").add_metabolites({Rubisco_balance:3})
    cobra_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1").add_metabolites({Rubisco_balance:-1})

    #generic ATPase and NADPH oxidase
    Maintenance_constraint = Metabolite("ATPase_NADPHoxidase_constraint_c1",name =  "ATPase_NADPHoxidase_constraint_c1", compartment = "c1")
    Maintenance_constraint2 = Metabolite("ATPase_NADPHoxidase_constraint_c2",name =  "ATPase_NADPHoxidase_constraint_c2", compartment = "c2")
    Maintenance_constraint3 = Metabolite("Light_dark_maintainence_constraint",name =  "Light_dark_maintainence_constraint", compartment = "c1")
    cobra_model.reactions.get_by_id("ATPase_tx1").add_metabolites({Maintenance_constraint:1,Maintenance_constraint3:1})
    cobra_model.reactions.get_by_id("ATPase_tx2").add_metabolites({Maintenance_constraint2:1,Maintenance_constraint3:-1})
    cobra_model.reactions.get_by_id("NADPHoxc_tx1").add_metabolites({Maintenance_constraint:-3})
    cobra_model.reactions.get_by_id("NADPHoxc_tx2").add_metabolites({Maintenance_constraint2:-3})
    cobra_model.reactions.get_by_id("NADPHoxm_tx1").add_metabolites({Maintenance_constraint:-3})
    cobra_model.reactions.get_by_id("NADPHoxm_tx2").add_metabolites({Maintenance_constraint2:-3})
    cobra_model.reactions.get_by_id("NADPHoxp_tx1").add_metabolites({Maintenance_constraint:-3})
    cobra_model.reactions.get_by_id("NADPHoxp_tx2").add_metabolites({Maintenance_constraint2:-3})

    ##constrain sucrose and starch storage
    Sucorse_starch_balance = Metabolite("sucrose_starch_bal_c", name = "Weights to balance sucrose-starch uptake", compartment = "c1")
    cobra_model.reactions.get_by_id("SUCROSE_v_dielTransfer").add_metabolites({Sucorse_starch_balance:-90})
    cobra_model.reactions.get_by_id("STARCH_p_dielTransfer").add_metabolites({Sucorse_starch_balance:10})

    #Plastid enolase was not detected in Arabidopsis mesophyll tissue
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p1").lower_bound=0
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p1").upper_bound=0
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p2").lower_bound=0
    cobra_model.reactions.get_by_id("2PGADEHYDRAT_RXN_p2").upper_bound=0

    #Setting chloroplastic NADPH dehydrogenase to 0  ((Yamamoto et al., 2011)
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").lower_bound=0
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p1").upper_bound=0
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").lower_bound=0
    cobra_model.reactions.get_by_id("NADPH_Dehydrogenase_p2").upper_bound=0

    #Set biomass to zero
    cobra_model.reactions.get_by_id("Biomass_tx1").lower_bound=0
    cobra_model.reactions.get_by_id("Biomass_tx1").upper_bound=0
    cobra_model.reactions.get_by_id("Biomass_tx2").lower_bound=0
    cobra_model.reactions.get_by_id("Biomass_tx2").upper_bound=0

    #ATP_ADP_Pi constrained to 0 because while there is evidence for its existance, it does not carry high flux
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc1").lower_bound = 0
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc1").upper_bound = 0
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc2").lower_bound = 0
    cobra_model.reactions.get_by_id("ATP_ADP_Pi_pc2").upper_bound = 0

    #turn off chlorophyll a/b cycling for energy dissipation
    cobra_model.reactions.get_by_id("RXN_7674_p1").lower_bound = 0
    cobra_model.reactions.get_by_id("RXN_7674_p1").upper_bound = 0

    #turn off cytosolic ethanol-ethanal cycle for NADH dissipation
    cobra_model.reactions.get_by_id("RXN_10745_NAD_c1").lower_bound = 0
    cobra_model.reactions.get_by_id("RXN_10745_NAD_c1").upper_bound = 0

    #turn off cytosolic ferric chelate reductase cycle for NADH dissipation
    cobra_model.reactions.get_by_id("FERRIC_CHELATE_REDUCTASE_RXN_c1").lower_bound = 0
    cobra_model.reactions.get_by_id("FERRIC_CHELATE_REDUCTASE_RXN_c1").upper_bound = 0

    #Adding a H_mc reaction to allow protons into mitochondria
    for i in range(1,3):
        rxn = Reaction("H_mc"+str(i))
        rxn.add_metabolites({cobra_model.metabolites.get_by_id("PROTON_c"+str(i)):-1,cobra_model.metabolites.get_by_id("PROTON_m"+str(i)):1})
        rxn.lower_bound=0
        rxn.upper_bound=1000
        cobra_model.add_reactions({rxn})

    return cobra_model

#####################################################################
# This function generates ATP budgets for a given flux distribution #
# inputs: 1) an FBA model, 2) a dictionary object with reaction ids #
# as keys and reaction fluxes as values, 3) name of output file (op-#
# -tional), 4) Option to show plots, 5) If choosing to show plot, c-#
# -hoose wether to use percentage or absolute values in the plot.   #
#####################################################################
def generateATPbudget(model,solution,outfile="",show_plot=True,percentage=False,day_or_night_tag="1",save_plot_to="temp.png"):
    if outfile!="":
        fout = open(outfile,"w")
    ATPdict = dict()
    total = 0
    for p in ("c","p","m","x"):
        met=model.metabolites.get_by_id("ATP_"+p+day_or_night_tag)
        met1=model.metabolites.get_by_id("aATP_"+p+day_or_night_tag)
        for rxn in met.reactions:
            if rxn.id.__contains__("ATP_AMP_mc") or rxn.id.__contains__("ATP_ADP_mc") or rxn.id.__contains__("ATP_pc") or rxn.id.__contains__("AMP_ATP_xc") or rxn.id.__contains__("ATP_ADP_Pi_pc"):
                continue
            sto=rxn.metabolites.get(met)
            sto1=rxn.metabolites.get(met1)
            if outfile!="":
                fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution.get(rxn.id)*(sto+sto1))+"\t"+met.compartment+"\n")
            ATPdict[rxn.id]=solution.get(rxn.id)*(sto+sto1)
            if solution.get(rxn.id)*(sto+sto1) > 0:
                total = total + (solution.get(rxn.id)*(sto+sto1))
    if outfile!="":
        fout.close()

    ATPdict2 = dict()
    ATPdict2["Others-pos"]=0
    ATPdict2["Others-neg"]=0
    baseline = dict()
    pos_base=0
    neg_base=0
    i=0
    for rxn in ATPdict.keys():
        if ATPdict[rxn]>0:
            if ATPdict[rxn] < total*0.05:
                if percentage:
                    ATPdict2["Others-pos"]=ATPdict2["Others-pos"]+float(ATPdict[rxn]*100)/total
                else:
                    ATPdict2["Others-pos"]=ATPdict2["Others-pos"]+ATPdict[rxn]
                continue
            base = pos_base
            if percentage:
                ATPdict2[rxn]=float(ATPdict[rxn]*100)/total
                pos_base = pos_base + float(ATPdict[rxn]*100)/total
            else:
                pos_base = pos_base + ATPdict[rxn]
                ATPdict2[rxn]=ATPdict[rxn]
        else:
            if abs(ATPdict[rxn]) < total*0.05:
                if percentage:
                    ATPdict2["Others-neg"]=ATPdict2["Others-neg"]+float(ATPdict[rxn]*100)/total
                else:
                    ATPdict2["Others-neg"]=ATPdict2["Others-neg"]+ATPdict[rxn]
                continue
            base = neg_base
            if percentage:
                ATPdict2[rxn]=float(ATPdict[rxn]*100)/total
                neg_base = neg_base + float(ATPdict[rxn]*100)/total
            else:
                neg_base = neg_base + ATPdict[rxn]
                ATPdict2[rxn]=ATPdict[rxn]
        i=i+1
        baseline[rxn]=base
    baseline["Others-pos"]=pos_base
    baseline["Others-neg"]=neg_base

    if show_plot:
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 10}) #sets a global fontsize
        plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
        plt.rcParams['xtick.major.width'] = 1
        plt.rcParams['ytick.major.size'] = 5
        plt.rcParams['ytick.major.width'] = 1
        plt.rcParams['axes.linewidth']=2 # makes axes line thicker
        plt.figure(figsize=(3,4))
        for rxn in ATPdict2.keys():
            plt.bar(1,ATPdict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn)
        plt.xlim(0.8,1.2)
        if percentage:
            plt.ylabel("ATP produced/consumed (%)")
        else:
            plt.ylabel("ATP produced/consumed (in moles)")
        handles, labels = plt.gca().get_legend_handles_labels()
        labels2=list(set(labels)-set(["Others-neg","Others-pos"]))+list(["Others-neg","Others-pos"])
        handles2=[handles[labels.index(i)] for i in labels2]
        lgd=plt.legend(handles2,labels2,bbox_to_anchor=(1,1))
        plt.axhline(0,linestyle="--",color="black")
        plt.tight_layout
        plt.savefig(save_plot_to, bbox_extra_artists=(lgd,), bbox_inches='tight')

#####################################################################
# This function generates ATP budgets for a given flux distribution #
# inputs: 1) an FBA model, 2) a dictionary object with reaction ids #
# as keys and reaction fluxes as values, 3) name of output file (op-#
# -tional), 4) Option to show plots, 5) If choosing to show plot, c-#
# -hoose wether to use percentage or absolute values in the plot 6) #
# Provide a day or night indicator tag to specify day or night NAD(-#
# -P)H summary 7) a destination file to save plot to                #
#####################################################################
def generateNADHNADPHbudget(model,solution,outfile="",show_plot=True,percentage=False,day_or_night_tag="1",save_plot_to="temp"):
    if outfile!="":
        fout = open(outfile,"w")
    Reddict = dict()
    total = 0
    for red in ["NADPH","NADH"]:
        for p in ("c","p","m","x"):
            if len(model.metabolites.query(red+"_"+p+day_or_night_tag))==0:
                continue
            met=model.metabolites.get_by_id(red+"_"+p+day_or_night_tag)
            for rxn in met.reactions:
                sto=rxn.metabolites.get(met)
                sto1=0#rxn.metabolites.get(met1)
                if outfile!="":
                    fout.write(rxn.id+"\t"+rxn.reaction+"\t"+str(solution.get(rxn.id)*(sto+sto1))+"\t"+met.compartment+"\n")
                Reddict[rxn.id]=solution.get(rxn.id)*(sto+sto1)
                if solution.get(rxn.id)*(sto+sto1) > 0:
                    total = total + (solution.get(rxn.id)*(sto+sto1))
    if outfile!="":
        fout.close()

    Reddict2 = dict()
    Reddict2["Others-pos"]=0
    Reddict2["Others-neg"]=0
    baseline = dict()
    pos_base=0
    neg_base=0
    i=0
    for rxn in Reddict.keys():
        if Reddict[rxn]>0:
            if Reddict[rxn] < total*0.05:
                if percentage:
                    Reddict2["Others-pos"]=Reddict2["Others-pos"]+float(Reddict[rxn]*100)/total
                else:
                    Reddict2["Others-pos"]=Reddict2["Others-pos"]+Reddict[rxn]
                continue
            base = pos_base
            if percentage:
                Reddict2[rxn]=float(Reddict[rxn]*100)/total
                pos_base = pos_base + float(Reddict[rxn]*100)/total
            else:
                pos_base = pos_base + Reddict[rxn]
                Reddict2[rxn]=Reddict[rxn]
        else:
            if abs(Reddict[rxn]) < total*0.05:
                if percentage:
                    Reddict2["Others-neg"]=Reddict2["Others-neg"]+float(Reddict[rxn]*100)/total
                else:
                    Reddict2["Others-neg"]=Reddict2["Others-neg"]+Reddict[rxn]
                continue
            base = neg_base
            if percentage:
                Reddict2[rxn]=float(Reddict[rxn]*100)/total
                neg_base = neg_base + float(Reddict[rxn]*100)/total
            else:
                neg_base = neg_base + Reddict[rxn]
                Reddict2[rxn]=Reddict[rxn]
        i=i+1
        baseline[rxn]=base
    baseline["Others-pos"]=pos_base
    baseline["Others-neg"]=neg_base

    if show_plot:
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size': 10}) #sets a global fontsize
        plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
        plt.rcParams['xtick.major.width'] = 1
        plt.rcParams['ytick.major.size'] = 5
        plt.rcParams['ytick.major.width'] = 1
        plt.rcParams['axes.linewidth']=2 # makes axes line thicker
        plt.figure(figsize=(3,4))
        for rxn in Reddict2.keys():
            plt.bar(1,Reddict2[rxn],width=0.1,bottom=baseline[rxn],label=rxn)
        plt.xlim(0.8,1.2)
        if percentage:
            plt.ylabel("NAD(P)H produced/consumed (%)")
        else:
            plt.ylabel("NAD(P)H produced/consumed (in moles)")
        handles, labels = plt.gca().get_legend_handles_labels()
        labels2=list(set(labels)-set(["Others-neg","Others-pos"]))+list(["Others-neg","Others-pos"])
        handles2=[handles[labels.index(i)] for i in labels2]
        lgd=plt.legend(handles2,labels2,bbox_to_anchor=(1,1))
        plt.axhline(0,linestyle="--",color="black")
        plt.tight_layout
        plt.savefig(save_plot_to, bbox_extra_artists=(lgd,), bbox_inches='tight')
