Aneurysm Simulation

This project uses the aneurysm tutorial dataset available at: http://www.visitusers.org/index.php?title=Tutorial_Data

VisIt software available at: https://wci.llnl.gov/simulation/computer-codes/visit/downloads


The aneurysm data set defines a mesh representing the affected blood vessel and a vector field showing the 3D velocity at each point. VisIt is easily able to plot the mesh and vector field, but I wanted to do something more. I had the idea to advect blood cells through the vessel as they are being pushed around by the velocity field. I used pathlines to see where a blood cell would travel given some initial position. The intermediate steps in these pathlines were saved to a file. This allowed me to simply iterate over each location and move each blood cell to its next position. 


instructions: 

	1. Modify hardcoded file paths in initialize.py as appropriate. 
	2. Launch Python command line interface from VisIt.  
	3. Source("path/to/initialize.py")
	4. execute main() to plot the data 


	
