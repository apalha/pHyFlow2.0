/**********************************************************************
* Description: NS ONLY meshing of the NACA0012
*	Gmsh Mesh for a airfoil
*
* References:
*	Gmsh Manual
*
* Copyright 2013. Lento Manickathan - l.manickathan@student.tudelft.nl
***********************************************************************/


Include "NACA0012.geo";



/******************************
*	Define Geometry
******************************/

// Length Scales
DefineConstant[ lcFar = { 2, Path "Gmsh/Parameters"}];

// Points - Far Field
pMid1 = newp; Point(pMid1) = {1, 0, 0, 0.5};
pMid2 = newp; Point(pMid2) = {2, 0, 0, 0.5};
pMid3 = newp; Point(pMid3) = {4, 0, 0, 0.5};



// Points - pExt
pMid1 = newp; Point(pMid1) = {4, 0, 0, 0.05};
pMid2 = newp; Point(pMid2) = {2, 1, 0, 0.05};
pMid3 = newp; Point(pMid3) = {-1, 0, 0, 0.05};
pMid4 = newp; Point(pMid4) = {2, -1, 0, 0.05};

pCenter2 = newp; Point(pCenter2) = {2, 0, 0, 0.05};

// Ellipse - lExt
lMid1 = newreg; Ellipse(lMid1) = {pMid1,pCenter2,pMid1,pMid2};
lMid2 = newreg; Ellipse(lMid2) = {pMid2,pCenter2,pMid3,pMid3};
lMid3 = newreg; Ellipse(lMid3) = {pMid3,pCenter2,pMid3,pMid4};
lMid4 = newreg; Ellipse(lMid4) = {pMid4,pCenter2,pMid1,pMid1};
lMidLL = newll; Line Loop(lMidLL) = {lMid1,lMid2,lMid3,lMid4};

// Points Far Field 
pFar1 = newp; Point(pFar1) = {10, 0, 0, lcFar};
pFar2 = newp; Point(pFar2) = {10, 8, 0, lcFar};
pFar3 = newp; Point(pFar3) = {-8, 8, 0, lcFar};
pFar4 = newp; Point(pFar4) = {-8, -8, 0, lcFar};
pFar5 = newp; Point(pFar5) = {10, -8, 0, lcFar};

lFar1 = newreg; Line(lFar1) = {pMid1, pFar1};
lFar2 = newreg; Line(lFar2) = {pFar1,pFar2};
lFar3 = newreg; Line(lFar3) = {pFar2,pFar3};
lFar4 = newreg; Line(lFar4) = {pFar3,pFar4};
lFar5 = newreg; Line(lFar5) = {pFar4,pFar5};
lFar6 = newreg; Line(lFar6) = {pFar5,pFar1};
lFarLL = newll; Line Loop(lFarLL) = {lFar1,lFar2,lFar3,lFar4,lFar5,lFar6};

// Transfinite Lines
Transfinite Line {-lMid2, lMid3} = 150 Using Progression 1.006;
Transfinite Line {lMid1, lMid4} = 50 Using Progression 1;
Transfinite Line {lFar1} = 50 Using Progression 1.03;
Transfinite Line {lFar2, -lFar6} = 25 Using Progression 1.05;
Transfinite Line {lFar4} = 20 Using Progression 1;

/******************************
*	Define Mesh
******************************/
fluidMid = newsl; Plane Surface(fluidMid) = {lExtLL, lMidLL};
fluidFar = newsl; Plane Surface(fluidFar) = {lFarLL,lMidLL};



// Mesh Parameters
meshAlgorithm = 8;  // delaunay quadrangulation
Mesh.Smoothing = 1;

// Physical Groups
Physical Surface(1) = {fluid, fluidFar, fluidMid};
Physical Line(2) = {lBody_top,lBody_bot};
Physical Line(3) = {lFar4,lFar3,lFar5};
Physical Line(4) = {lFar2,lFar6};

