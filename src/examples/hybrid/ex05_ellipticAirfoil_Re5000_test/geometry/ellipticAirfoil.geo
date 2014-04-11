/**********************************************************************
* Description: 
*	Gmsh Mesh for a elliptic airfoil
*
* References:
*	Gmsh Manual
*
* Copyright 2013. Lento Manickathan - l.manickathan@student.tudelft.nl
***********************************************************************/

/******************************
*	Define Geometry
******************************/

// Length Scales
DefineConstant[ lcBody  = { 0.00333333333333, Path "Gmsh/Parameters"}];
DefineConstant[ lcExt  = { 0.00666666666667, Path "Gmsh/Parameters"}];

xCM = 0.25;

// Points - pCenter
pCenter1 = newp; Point(pCenter1) = {0 + xCM, 0, 0, 0.00333333};



// Points - pBody
pBody1 = newp; Point(pBody1) = {0.5 + xCM, 0, 0, 0.00333333333333};
pBody2 = newp; Point(pBody2) = {0 + xCM, 0.05, 0, 0.00333333333333};
pBody3 = newp; Point(pBody3) = {-0.5 + xCM, 0, 0, 0.00333333333333};
pBody4 = newp; Point(pBody4) = {0 + xCM, -0.05, 0, 0.00333333333333};

// Points - pExt
pExt1 = newp; Point(pExt1) = {0.625 + xCM, 0, 0, 0.00666667};
pExt2 = newp; Point(pExt2) = {0 + xCM, 0.375, 0, 0.00666667};
pExt3 = newp; Point(pExt3) = {-0.625 + 0.25, 0, 0, 0.00666667};
pExt4 = newp; Point(pExt4) = {0 + xCM, -0.375, 0, 0.00666667};




// Ellipse - lExt
lBody1 = newreg; Ellipse(lBody1) = {pBody1,pCenter1,pBody1,pBody2};
lBody2 = newreg; Ellipse(lBody2) = {pBody2,pCenter1,pBody3,pBody3};
lBody3 = newreg; Ellipse(lBody3) = {pBody3,pCenter1,pBody3,pBody4};
lBody4 = newreg; Ellipse(lBody4) = {pBody4,pCenter1,pBody1,pBody1};
lBodyLL = newll; Line Loop(lBodyLL) = {lBody1,lBody2,lBody3,lBody4};

// Ellipse - lExt
lExt1 = newreg; Ellipse(lExt1) = {pExt1,pCenter1,pExt1,pExt2};
lExt2 = newreg; Ellipse(lExt2) = {pExt2,pCenter1,pExt3,pExt3};
lExt3 = newreg; Ellipse(lExt3) = {pExt3,pCenter1,pExt3,pExt4};
lExt4 = newreg; Ellipse(lExt4) = {pExt4,pCenter1,pExt1,pExt1};
lExtLL = newll; Line Loop(lExtLL) = {lExt1,lExt2,lExt3,lExt4};


Transfinite Line {lExt1,-lExt2,lExt3,-lExt4} = 150 Using Progression 1;
Transfinite Line {lBody1,-lBody2,lBody3,-lBody4} = 200 Using Progression 1.003;
//Transfinite Line {lExt1,-lExt2,lExt3,-lExt4} = 150 Using Progression 1;
/******************************
*	Define Mesh
******************************/

// Mesh Domains
fluid = newsl; Plane Surface(fluid) = {lBodyLL,lExtLL};
Mesh.Smoothing = 100;

// Physical Groups
Physical Surface(1) = {fluid};
Physical Line(2)    = {lBody1,lBody2,lBody3,lBody4};
Physical Line(3)    = {lExt1,lExt2,lExt3,lExt4};

