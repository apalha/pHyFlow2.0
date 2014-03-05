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
DefineConstant[ lcFar   = { 5.0, Path "Gmsh/Parameters"}]; 		// Outer Domain
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
**** Defining the outer region
************************************/


// Outer Domain parameters
inner_xmin = x0-4*R; inner_xmax = x0+20*R;
inner_ymin = y0-4*R; inner_ymax = y0+4*R;

inner2_xmin = x0-12*R; inner2_xmax = x0+60*R;
inner2_ymin = y0-12*R; inner2_ymax = y0+12*R;



// Points
inner_p1 = newp; Point(inner_p1) = {inner_xmin, inner_ymin, 0, lcCoarse*5};
inner_p2 = newp; Point(inner_p2) = {inner_xmin, inner_ymax, 0, lcCoarse*5};
inner_p3 = newp; Point(inner_p3) = {inner_xmax, inner_ymax, 0, lcCoarse*5};
inner_p4 = newp; Point(inner_p4) = {inner_xmax, inner_ymin, 0, lcCoarse*5};


// Points
inner2_p1 = newp; Point(inner2_p1) = {inner2_xmin, inner2_ymin, 0, lcCoarse*30};
inner2_p2 = newp; Point(inner2_p2) = {inner2_xmin, inner2_ymax, 0, lcCoarse*30};
inner2_p3 = newp; Point(inner2_p3) = {inner2_xmax, inner2_ymax, 0, lcCoarse*30};
inner2_p4 = newp; Point(inner2_p4) = {inner2_xmax, inner2_ymin, 0, lcCoarse*30};

// Lines
//inner_L12 = newl; Line(inner_L12) = {inner_p1, inner_p2};
//inner_L23 = newl; Line(inner_L23) = {inner_p2, inner_p3};
//inner_L34 = newl; Line(inner_L34) = {inner_p3, inner_p4};
//inner_L41 = newl; Line(inner_L41) = {inner_p4, inner_p1};

//Transfinite Line {inner_L23,-inner_L41} = 150 Using Progression 1.006;
//Transfinite Line {inner_L34} = 30 Using Progression 1;
//inner_Box = newll; Line Loop(inner_Box) = {inner_L12, inner_L23, inner_L34, inner_L41};

splineLine = newreg; Spline(splineLine) = {inner_p3, inner_p2,inner_p1, inner_p4, inner_p3};
splineLine2 = newreg; Spline(splineLine2) = {inner2_p3, inner2_p2,inner2_p1, inner2_p4, inner2_p3};
//Transfinite Line {splineLine} = 400 Using Bump 0.8;

inner_Box = newll; Line Loop(inner_Box) = {splineLine};
inner2_Box = newll; Line Loop(inner2_Box) = {splineLine2};


// Outer Domain parameters
outer_xmin = x0-60*R; outer_xmax = x0+100*R;
outer_ymin = y0-60*R; outer_ymax = y0+60*R;

// Points
outer_p1 = newp; Point(outer_p1) = {outer_xmin, outer_ymin, 0, lcFar};
outer_p2 = newp; Point(outer_p2) = {outer_xmin, outer_ymax, 0, lcFar};
outer_p3 = newp; Point(outer_p3) = {outer_xmax, outer_ymax, 0, lcFar};
outer_p4 = newp; Point(outer_p4) = {outer_xmax, outer_ymin, 0, lcFar};

// Lines
outer_L12 = newl; Line(outer_L12) = {outer_p1, outer_p2};
outer_L23 = newl; Line(outer_L23) = {outer_p2, outer_p3};
outer_L34 = newl; Line(outer_L34) = {outer_p3, outer_p4};
outer_L41 = newl; Line(outer_L41) = {outer_p4, outer_p1};

Transfinite Line {outer_L34, outer_L12} = 50 Using Bump 2;

outer_Box = newll; Line Loop(outer_Box) = {outer_L12, outer_L23, outer_L34, outer_L41};


/************************************
**** Mesh Domains
************************************/
ext_Domain = newsl; Plane Surface(ext_Domain) = {body_Box, ext_Box};
inner_Domain = newsl; Plane Surface(inner_Domain) = {inner_Box, ext_Box};
inner2_Domain = newsl; Plane Surface(inner2_Domain) = {inner2_Box, inner_Box};
outer_Domain = newsl; Plane Surface(outer_Domain) = {outer_Box, inner2_Box};

Mesh.Algorithm = 8;
Mesh.Smoothing = 100;

/************************************
**** Physical Groups - Externel references
************************************/

// Physical Groups - Naming the boundaries
Physical Surface(1) = {ext_Domain, inner_Domain, inner2_Domain, outer_Domain};         // Fluid Domain
Physical Line(2)    = {L12, L21}; 		   // body boundary [2]
Physical Line(3) = {outer_L12,outer_L23, outer_L41};						// in-flow [3]
Physical Line(4) = {outer_L34};						// out-flow [4]
