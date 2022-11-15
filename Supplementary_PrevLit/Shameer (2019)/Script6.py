#######Observe mATP synthesis in growing C3 leaf#######
#In this script, a fully expanded C3 source leaf with #
#plastidic ATP shuttles constrained based on published#
#studies is modelled in a scenario when photon uptake #
#not forced on the system. Fluxes are observed for Fi-#
#gure 5c                                              #
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

#constrain chloroplast ATP shuttles based on published data
#For PEP-Pyr shuttle, constrain cplast PPDK based on Ishimaru, K., Ichikawa, H., Matsuoka, M., & Ohsugi, R. (1997). Analysis of a C4 maize pyruvate, orthophosphate dikinase expressed in C3 transgenic Arabidopsis plants. Plant Science, 129(1), 57-64.
cobra_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p1").lower_bound = 0
cobra_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p1").upper_bound = 0.03375
#For TP-3PGA shuttle, constrain cytosolic phosphorylating GAPDH activity based on Gibon, Y., et al(2004). A robot-based platform to measure multiple enzyme activities in Arabidopsis using a set of cycling assays: comparison of changes of enzyme activities and transcript levels during diurnal cycles and in prolonged darkness. The Plant Cell, 16(12), 3304-3325.
cobra_model.reactions.get_by_id("GAPOXNPHOSPHN_RXN_c1").lower_bound = -95.387597
cobra_model.reactions.get_by_id("GAPOXNPHOSPHN_RXN_c1").upper_bound = 95.387597
#and constrain cytosolic non-phosphorylating GAPDH based on Strand, Å., Zrenner, R., Trevanion, S., Stitt, M., Gustafsson, P., & Gardeström, P. (2000). Decreased expression of two key enzymes in the sucrose biosynthesis pathway, cytosolic fructose‐1, 6‐bisphosphatase and sucrose phosphate synthase, has remarkably different consequences for photosynthetic carbon metabolism in transgenic Arabidopsis thaliana. The Plant Journal, 23(6), 759-770.
cobra_model.reactions.get_by_id("1_PERIOD_2_PERIOD_1_PERIOD_9_RXN_c1").lower_bound = 0
cobra_model.reactions.get_by_id("1_PERIOD_2_PERIOD_1_PERIOD_9_RXN_c1").upper_bound = 0.33

#create a backup of source C3 model
phloem_model = cobra_model.copy()

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

#Model flux distribution at 200 PPFD in based on A vs PPFD data
WTphloemOut = dict()
solutiondict_phloem = dict()
FVAdict_phloem = dict()
PPFD_list=list()

PPFD = 200
tempModel2 = phloem_model.copy()
# constrain photon uptake upper bound
tempModel2.reactions.get_by_id("Photon_tx1").upper_bound = PPFD*0.9
#tempModel2.reactions.get_by_id("Photon_tx1").lower_bound = PPFD*0.9
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
PPFD_list.append(PPFD)
#perform FVA
fva_result = flux_analysis.flux_variability_analysis(tempModel2,fraction_of_optimum=1)
#record FVA
FVAdict_phloem[PPFD] = fva_result

#Print pFBA and FVA results
printFluxesAndFVA(tempModel2,solutiondict_phloem,FVAdict_phloem,outfile="Fluxes_Source_Leaf_ATPshuttlesconstrained_LightNotForced.csv")
