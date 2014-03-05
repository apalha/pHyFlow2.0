/**********************************************************************
* Description: 
*	Gmsh Mesh for a [0,1] x [0,1] 2D rectangle
* 	 - Gmsh project created on Wed Aug  15
*
* References:
*	Gmsh Manual
*
* Copyright 2013. Lento Manickathan - l.manickathan@student.tudelft.nl
***********************************************************************/


/************************************
**** Mesh parameteres, constants 
************************************/

// Length Scales
dX = 1; dY = 0.5;
x0 = 0.0; y0 = 0.0;

DefineConstant[ lc_X = {0.02, Path "Gmsh/Parameters"}]; //dX/100
DefineConstant[ lc_Y = {0.02, Path "Gmsh/Parameters"}]; //dY/50


// Domain definition
xmin = x0 - dX*0.5; xmax = x0 + dX*0.5;
ymin = y0 - dY*0.5; ymax = y0 + dY*0.5;


/************************************
**** Inner
************************************/

// Points
p1 = newp; Point(p1) = {xmin,   ymin, 0, lc_X};
p2 = newp; Point(p2) = {xmin,   ymax, 0, lc_X};
p3 = newp; Point(p3) = {0,      ymax, 0, lc_X};
p4 = newp; Point(p4) = {xmax,   ymax, 0, lc_X};
p5 = newp; Point(p5) = {xmax,   ymin, 0, lc_X};
p6 = newp; Point(p6) = {0,      ymin, 0, lc_X};


// Lines
L12 = newl; Line(L12) = {p1, p2};
L23 = newl; Line(L23) = {p2, p3};
L36 = newl; Line(L36) = {p3, p6};
L61 = newl; Line(L61) = {p6, p1};

L34 = newl; Line(L34) = {p3, p4};
L45 = newl; Line(L45) = {p4, p5};
L56 = newl; Line(L56) = {p5, p6};



Box_inner = newll; Line Loop(Box_inner) = {L12, L23, L36, L61};
Box_inner2 = newll; Line Loop(Box_inner2) = {-L36, L34, L45, L56};

//Transfinite Line {L36} = 50 Using Progression 1;
Transfinite Line {L12,L23,L36,L61,L34,L45,L56} = 50 Using Progression 1;

/************************************
**** Mesh Domains
************************************/

Domain_inner = newsl; Plane Surface(Domain_inner) = {Box_inner};
Domain_inner2 = newsl; Plane Surface(Domain_inner2) = {Box_inner2};

Mesh.Smoothing = 100;

/************************************
**** Physical Groups - Externel references
************************************/


// Physical Groups - Naming the boundaries
Physical Surface(1) = {Domain_inner, Domain_inner2}; // Fluid Domain
Physical Line(3) = {L12, L23,  L34, L45, L56, L61};  // exterior boundary [3]

