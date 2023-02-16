# Shameer-et-al-Role-of-mitochondria-in-C3-leaf-during-the-day
This repository contains scripts used to support our study the role of mitochondria in C3 leaves during the light-phase

## Dependencies:

libsbml version 5.17.0  
pandas version 0.23.4  
[cobrapy version 0.13.4](https://github.com/opencobra/cobrapy/releases/tag/0.13.4)

## Instructions:
1. Install python version 2, libsbml, pandas and cobrapy
2. Go to folder contaning scripts
3. Open a console(max), terminal(linux) or command prompt(windows)
4. Run sripts (tw options available)
..a) using the command `python scriptname.py` where scriptname.py is the filename of the script
<p align='center'>
  OR
</p>
<p>
  &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp b) Enter python terminal using the command "python" and copy paste script contents to the python terminal
</p>

  
The latter method is recommended as it maintains the data in the python environment allowing one to analyse metabolic fluxes directly and make your own figures.

## Index:
* Script1.py - estiamting curve used to calculate CO2 assimilation rate from photosynthetic photon flux density(PPFD) based on data from Donahue et al 1997  
* Script1.ipynb - a jupyter notebok version of Script1.py   
* Script2.py - generating Figure 3a  
* Script2.ipynb - a jupyter notebok version of Script2.py   
* Script3.py - generating Figure 3b  
* Script3.ipynb - a jupyter notebok version of Script3.py   
* Script4.py - generating data for Figure 4  
* Script4.ipynb - a jupyter notebok version of Script4.py   
* Script5.py - generating Figure 5a and 5b  
* Script5.ipynb - a jupyter notebok version of Script5.py   
* Script6.py - generating data for Figure 5c  
* Script6.ipynb - a jupyter notebok version of Script6.py   
* Script7.py - generating data for Figure 5d  
* Script7.ipynb - a jupyter notebok version of Script7.py   
* studyFunctions - a python module with functions used in scripts 1-7   
* core_model.xml - a model representing primary metabolism in plant cells  
* Gomes2015_Biomass.csv - a text file with Arabidopsis biomass composition gathered from Gomes et al 2015  
* MetabolitesToTransfer.txt - a text file listing metabolites that will be allowed to accumulate over the diel cycle  
