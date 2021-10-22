To run, you can type into a terminal:


python fillgaps.py("filename.csv","results.csv")

where 
"filename.csv" is the file containing the data with missing points formatted as comma-delimited two-column csv
               in the following format:
	       <time data, numeric, equispaced>, <data with missing points>
"results.csv" is the name of the file you wish to save the imputed results to






You can also import the module in a python environment

import fillgaps as fg
fg.fillgaps("filename.csv","results.csv")