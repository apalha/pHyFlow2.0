/**********************************************************************
* Description: 
*	Gmsh Mesh for a 2D cylinder
* 	 - Gmsh project created on Wed Aug  25 13:37:25 2013
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
DefineConstant[ lcCoarse= { 0.03, Path "Gmsh/Parameters"}]; 		// Outer Domain
DefineConstant[ lcFine  = { 0.01, Path "Gmsh/Parameters"}];		// Inner Refined Domain

// Body definition
x0 = 0.0; 
y0 = 0.0;
R  = 1.0;

pc = newp; Point(pc) = {0, 0, 0, lcFine};

// Cylinder Parameters
xmin = x0-R; xmax = x0+R;

// Exterior region parameters
ext = R*0.5;
xminExt = x0-(R+ext);
xmaxExt = x0+(R+ext);


/************************************
**** Body : body description 
************************************/

// Points
p1 = newp; Point(p1) = {xmin, 0,   0, lcFine};
p2 = newp; Point(p2) = {xmax, 0,   0, lcFine};

// Lines
L12 = newreg; Circle(L12) = {p1, pc, p2};
L21 = newreg; Circle(L21) = {p2, pc, p1};


//Transfinite Line {L12,L21} = 250 Using Bump 0.25;

body_Box = newll; Line Loop(body_Box) = {L12, L21};

/************************************
**** Defining the exterior region
************************************/

// Points
p1Ext = newp; Point(p1Ext) = {xminExt, 0, 0, lcCoarse};
p2Ext = newp; Point(p2Ext) = {xmaxExt, 0, 0, lcCoarse};

// Lines
L12Ext = newreg; Circle(L12Ext) = {p1Ext, pc, p2Ext};
L21Ext = newreg; Circle(L21Ext) = {p2Ext, pc, p1Ext};


// Line loop
ext_Box = newll; Line Loop(ext_Box) = {L12Ext, L21Ext};


/************************************
**** Mesh Domains
************************************/
ext_Domain = newsl; Plane Surface(ext_Domain) = {body_Box, ext_Box};

Mesh.Algorithm = 8;
Mesh.Smoothing = 100;

/************************************
**** Physical Groups - Externel references
************************finf************/

// Physical Groups - Naming the boundaries
Physical Surface(1) = {ext_Domain};         // Fluid Domain
Physical Line(2)    = {L12, L21}; 		   // body boundary [2]
Physical Line(3)    = {L12Ext, L21Ext};    // exterior-boundary

