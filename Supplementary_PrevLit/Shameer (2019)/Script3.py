#######Observe mATP synthesis in growing C3 leaf#######
#In this script, a growing C3 leaf is modelled and fl-#
#-ux through mitochondrial ATP synthase is observed to#
# generate Figure 3b                                  #
#######################################################

#import functions from library
from libsbml import readSBML
from cobra import io,flux_analysis
from cobra.core import Metabolite, Reaction
from studyFunctions import *

#import model. Update file name and location in the next line
cobra_model = io.sbml.create_cobra_model_from_sbml_file("core_model.xml")

#set up a diel model, allow for day-night metabolite accumulations and constrain model to represent C3 leaf diel metabolism
cobra_model=setupDielModel(cobra_model,"MetabolitesToTransfer.txt")

#increasing upper bounds of model for higher light intensity simulations (>1000)
for rxn in cobra_model.reactions:
    if rxn.lower_bound < -999:
        rxn.lower_bound = -2000
    if rxn.upper_bound > 999:
        rxn.upper_bound = 2000

#transform source C3 model to a growing C3 leaf model
Biomass_model = cobra_model.copy()
Biomass_model = updateBiomassComposition(Biomass_model,"Gomes2015_Biomass.csv")
for i in range(1,3):
    met = Biomass_model.metabolites.get_by_id("X_Phloem_contribution_t"+str(i))
    Biomass_model.reactions.get_by_id("Phloem_output_tx"+str(i)).add_metabolites({met:-1})
    Biomass_model.reactions.get_by_id("Phloem_output_tx"+str(i)).lower_bound = 0
    Biomass_model.reactions.get_by_id("Phloem_output_tx"+str(i)).upper_bound = 0
    Biomass_model.reactions.get_by_id("Biomass_tx"+str(i)).add_metabolites({met:1})
    Biomass_model.reactions.get_by_id("Biomass_tx"+str(i)).lower_bound = 0
    Biomass_model.reactions.get_by_id("Biomass_tx"+str(i)).upper_bound = 1000
    Biomass_model.metabolites.get_by_id("X_Phloem_contribution_t"+str(i)).id = "X_Biomass_contribution_t"+str(i)
    Biomass_model.metabolites.get_by_id("X_Biomass_contribution_t"+str(i)).name = "X_Biomass_contribution[t]"

Biomass_model.metabolites._generate_index()
Biomass_model.reactions.get_by_id("diel_biomass").name = "Diel biomass accumulation"


#Arabidopsis assimilation rate (A) vs photosynthetic photon flux density (PPFD)
#data was gathered from Donahue et al 1997, see Script1.py
#if x = light, y = net CO2 uptake, y = a + bx + c^2
a = 0.0871351015562
b = 0.0291441670197
c = -0.000009013134

light=list()
CO2=list()
for x in range(100,1550,50):
    light.append(x)
    y=((a)+(b*x)+(c*(x**2)))
    CO2.append(y)

light_CO2=dict(zip(light, CO2))
import pandas as pd
df = pd.DataFrame(data={"PPFD":light,"Net CO2 uptake":CO2})
df = df[["PPFD","Net CO2 uptake"]]

#Model flux distribution in based on A vs PPFD data
WTphloemOut = dict()
solutiondict_phloem = dict()
FVAdict_phloem = dict()
PPFD_list=list()
for i in range(100,1550,50):
    PPFD = i
    tempModel2 = Biomass_model.copy()
    # constrain photon uptake
    tempModel2.reactions.get_by_id("Photon_tx1").upper_bound = PPFD*0.9
    tempModel2.reactions.get_by_id("Photon_tx1").lower_bound = PPFD*0.9
    # constrain maintenance cost
    tempModel2.reactions.get_by_id("ATPase_tx1").upper_bound = estimateMaintenance(200)
    tempModel2.reactions.get_by_id("ATPase_tx1").lower_bound = estimateMaintenance(200)
    #check for fesible solution
    solution=flux_analysis.parsimonious.pfba(tempModel2)
    #identify WT diel phloem export rate
    WTphloemOut[PPFD] = estimateOutputFromNetCO2(tempModel2,light_CO2[PPFD],verbose=False)
    #constrain diel phloem export rate
    tempModel2.reactions.get_by_id("diel_biomass").upper_bound = round(WTphloemOut[PPFD],3)
    tempModel2.reactions.get_by_id("diel_biomass").lower_bound = round(WTphloemOut[PPFD],3)
    #perform pFBA
    solution=flux_analysis.parsimonious.pfba(tempModel2)
    #record flux distribution
    solutiondict_phloem[PPFD]=solution.x_dict
    #record PPFD
    PPFD_list.append(i)

#Plot mitochondrial ATP synthase and diel phloem export rate for various PPFD
xlist = list()
ylist1 = list()
ylist2 = list()
for PPFD in PPFD_list:
    xlist.append(PPFD)
    ylist1.append(solutiondict_phloem[PPFD]["Mitochondrial_ATP_Synthase_m1"])
    ylist2.append(WTphloemOut[PPFD]*4/1000*60*60*24)

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14}) #sets a global fontsize
plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['axes.linewidth']=2 # makes axes line thicker
fig, ax1 = plt.subplots(figsize=(4, 4))
ax1.plot(xlist,ylist1,"go-")
ax1.set_xlabel("PPFD ($\mu$mol m$^2$ s$^{-1}$)")
ax1.set_ylabel("Daytime mitochondrial ATP\nsynthase flux ($\mu$mol m$^2$ s$^{-1}$)",color="g")
ax1.tick_params("y",colors="g")
ax2 = ax1.twinx()
ax2.plot(xlist,ylist2,"bs-")
ax2.set_ylabel("phloem output (mmol m$^2$ 24h$^{-1}$)",color="b",rotation='270')
ax2.yaxis.set_label_coords(1.27,0.5)
ax2.tick_params("y",colors="b")
fig.tight_layout()
plt.show()

#print pFBA and FVA results to file
WTphloemOut = dict()
solutiondict_phloem = dict()
FVAdict_phloem = dict()
PPFD_list=list()

for i in range(100,1550,50):
    PPFD = i
    print("Running PPFD = "+str(i))
    tempModel2 = phloem_model.copy()
    tempModel2.reactions.get_by_id("Photon_tx1").upper_bound = PPFD*0.9
    tempModel2.reactions.get_by_id("Photon_tx1").lower_bound = PPFD*0.9
    tempModel2.reactions.get_by_id("ATPase_tx1").upper_bound = estimateMaintenance(200)
    tempModel2.reactions.get_by_id("ATPase_tx1").lower_bound = estimateMaintenance(200)
    solution=flux_analysis.parsimonious.pfba(tempModel2)
    WTphloemOut[PPFD] = estimateOutputFromNetCO2(tempModel2,light_CO2[PPFD],verbose=False)
    tempModel2.reactions.get_by_id("diel_biomass").upper_bound = round(WTphloemOut[PPFD],3)
    tempModel2.reactions.get_by_id("diel_biomass").lower_bound = round(WTphloemOut[PPFD],3)
    solution=flux_analysis.parsimonious.pfba(tempModel2)
    solutiondict_phloem[PPFD]=solution.x_dict
    PPFD_list.append(i)
    fva_result = flux_analysis.flux_variability_analysis(tempModel2,fraction_of_optimum=1)
    FVAdict_phloem[PPFD] = fva_result


printFluxesAndFVA(tempModel2,solutiondict_phloem,FVAdict_phloem,outfile="Fluxes_Growing_Leaf.csv")
print("Done")
