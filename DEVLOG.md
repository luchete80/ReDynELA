20241108 - Fixed remesh commands (was meshing already) 
         - Changed indices to begin with 1 like in example
         - Added gmsh export (like in example )
20241111 - Added Plastic Strain remesh criteria
         - Added ExplicitSolverRmsh class
20241122 - Copy fields to new Structure
20241125 - Fixed several things regarding to trias. 
         - ElementData. Shape Functions, internal gauss coordinates. 
         - Still remaining: 
         - element->computeGlob2Loc();
         - element->computeBoundBox();
--------------------------------------------------------
20241204 - ADDED tethex lib (as submodule, different from mmg which is external)
         - FIXED mmg (removed link to FE library)
         - See that in 5.8.0 version built with mingw has tochange some things)