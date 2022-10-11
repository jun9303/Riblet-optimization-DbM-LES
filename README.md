# Riblet-optimization-DbM-LES-BO
 
How to run:
1. compile the LES solver (cd solver && make new && cd ..). A Fortran compiler is required (recommend: Intel Fortran Compiler)
2. create the testspace (mkdir testspace)
3. run ./eval.py. 13 command line arguments are needed: w_1 w_2 ... w_10 w+ h+ s+. (./eval.py w_1 w_2 w_3 w_4 w_5 w_6 w_7 w_8 w_9 w_10 w+ h+ s+)
4. after the simulation ends, fribs.dat in the data folder contains the drag data on the channel walls
	- Column 1: time
	- Column 2: instantaneous drag on the riblet wall
	- Column 3: instantaneous drag on the flat wall
	- Column 4: time-averaged drag on the riblet wall
	- Column 5: time-averaged drag on the flat wall