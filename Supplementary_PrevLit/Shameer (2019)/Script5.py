#######Observe mATP synthesis in growing C3 leaf#######
#In this script, a fully expanded C3 source leaf with #
#plastidic ATP shuttles constrained based on published#
#studies is modelled in a scenario when photon uptake #
#not forced on the system. Fluxes are observed for Fi-#
#gure 5a and 5b                                       #
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

#Model flux distribution over a range of light intensities based on A vs PPFD data
WTphloemOut = dict()
solutiondict_phloem = dict()
FVAdict_phloem = dict()
PPFD_list=list()
mitoATPase=dict()
light_used=dict()

for i in range(100,1550,100):
    PPFD = i
    print("Running PPFD = "+str(i))
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
    PPFD_list.append(i)
    #record mATP synthesis rate
    mitoATPase[i]=tempModel2.reactions.get_by_id("Mitochondrial_ATP_Synthase_m1").x
    #record photon usage
    light_used[i]=tempModel2.reactions.get_by_id("Photon_tx1").x


#Generate plots for Figure 5a and b
import matplotlib.pyplot as plt
get_ipython().magic(u'matplotlib inline')
plt.rcParams.update({'font.size': 14}) #sets a global fontsize
plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['ytick.major.width'] = 2
plt.rcParams['axes.linewidth']=2 # makes axes line thicker

import pandas as pd
import numpy as np7

df=pd.DataFrame(data={'PPFD':PPFD_list})
df['Mito ATPase flux']= df['PPFD'].map(mitoATPase)
df['Phloem output flux'] =df['PPFD'].map(WTphloemOut)
df['Diel phloem output mmol/24h'] = df['Phloem output flux']*4/1000*60*60*24
df['CO2_tx1']=df['PPFD'].map(light_CO2)
df['light used']=df['PPFD'].map(light_used)
df['% light used']=df['light used']/df['PPFD']*100


fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8, 4), sharex=True)
df.plot(ax=ax[0], x='PPFD', y='% light used', legend=False, lw=2, linestyle='-', marker='s', color='black', label = 'Light used')
df.plot(ax=ax[1], x='PPFD', y='Mito ATPase flux', legend=False, lw=2, linestyle='-', marker='s', color='g',label ='M ATPase')
ax[0].set_xlabel('PPFD ($\mu$mol m$^2$ s$^{-1}$)',fontsize=12)
ax[0].set_ylabel('% of photons used for photochemistry',fontsize=12)
ax[1].set_xlabel('PPFD ($\mu$mol m$^2$ s$^{-1}$)',fontsize=12)
ax[1].set_ylabel('Daytime mitochondrial ATP\nsynthase flux ($\mu$mol m$^2$ s$^{-1}$)',fontsize=12)
ax[0].set_ylim(20, 50)
ax[1].set_xlim(0, 1550)
ax[1].set_ylim(0, 4.5)

#ax1.plot(np.nan, 'bs-', label = 'light used')  # Add a blank dataset to ax with correct line for ax2 legend
#ax1.legend(loc = 'lower right', bbox_to_anchor=(1.0, 0.05))
#plt.savefig('8_2_phloem_light_vs_mitoATPase_lightupperBonly.pdf', bbox_inches='tight')
plt.tight_layout()
plt.show()
