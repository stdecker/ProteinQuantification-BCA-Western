# ProteinQuantification-BCA-Western
Protein Quantification Following Standard BCA Assay. This app allows for quick quantification of protein concentration based on a standard BCA assay ([like these](https://www.thermofisher.com/us/en/home/life-science/protein-biology/protein-assays-analysis/protein-assays/bca-protein-assays.html)). Please note that this assay must be laid out as in the example below.

![image](https://github.com/stdecker/ProteinQuantification-BCA-Western/assets/35308658/03f31cf6-c395-4be5-b1e3-45897594a98e)

## Volume of Sample for Desired Amount of Protein

This assay is run using R/RStudio and Shiny. Please make sure all necessary packages are installed (readxl, shiny) using `install.packages('readxl', 'shiny')` before running the program. Once installed, highlight all (Ctrl + A) and run (Ctrl + Enter). Shiny should pop up. 

![image](https://github.com/stdecker/ProteinQuantification-BCA-Western/assets/35308658/8e43c936-bcb6-40dd-a05b-a15daa85b2e7)

From there, upload your BCA Assay results from Excel using the Browse button. Update the minimum and maximum standard concentrations, the number of standards used for the standard curve (number of rows used), and the difference in concentration between each step. For example, if you were to set up 6 wells of standards, ranging from 0 &mu;g/mL to 5 &mu;g/mL in 1 &mu;g/mL increments, you would put the following:
  Minumum: 0
  Maximum: 5
  Number of Standards: 6
  Steps: 1

The output will give you the volume of sample needed for the desired amount of protein per aliquot (default is 75 &mu;g, but can also be changed).

## Western Blot Mixture

This assay also allows for the quick construction of Western Blot Assay prep. To begin, follow the steps above. Change the amount of protein desired for each WB well, and the total volume desired for the WB well. If using RIPA, check the box and change volume of RIPA. Last, select the Laemmmli Buffer Concentration.
