/*********************************************************************
 * Ellipse perturbed by one corner on the major axis.
*
 *  - The Euler region is on the disk of radius R, but the PML is only on
 * the disk of radius R_PML with R_PML<R.
 *  - Unstructured mesh, structured at the sign-changing interface.
 *  - Corner at (x_c,0), symmetric w.r.t. x-axis
 *  - Corner disk discretized in Euler coordinates (z,theta)
 *  - Recombination algorithm can be used to produce quad mesh.
 */

/* Possible improvements:
 *   - Enforce mesh symmety w.r.t x-axis. Gmsh cannot enforce unstructured mesh symmetry,
 *  but adding a horizontal edge y=0 does the trick.
 *   - Implemenent mesh coarsening using background fields.
 */

SetFactory("OpenCASCADE");
    // ---- Input parameters

    // Ellipse gamma-m semi-axes
If (!Exists(a_m)) a_m=2.5; EndIf
If (!Exists(b_m)) b_m=1; EndIf
    // Ellipse gamma-d semi-axes
If (!Exists(a_d)) a_d=1.3*a_m; EndIf
If (!Exists(b_d)) b_d=Sqrt(a_d^2-Abs(a_m^2-b_m^2)); EndIf
    // Corner angle (rad)
If (!Exists(phi)) phi=Pi/2; EndIf
If (!Exists(R_PML)) R_PML = 0.8; EndIf // R_PML/R radius of PML region
If (!Exists(R_TR)) R_TR = 0.05; EndIf // R_TR/R radius of truncated region
    // Corner region offset
If (!Exists(x_ofst)) x_ofst = -a_d; EndIf
If (!Exists(y_ofst)) y_ofst = -b_d-Pi; EndIf
    // Mesh parameters (Structured interface)
    // outer ellipse
If (!Exists(a_o)) a_o=(1-0.1)*a_m+0.1*a_d; EndIf
If (!Exists(b_o)) b_o=Sqrt(Abs(a_o^2-Abs(a_m^2-b_m^2))); EndIf
    // inner ellipse
If (!Exists(a_i)) a_i=(1-2e-2)*a_m; EndIf    
If (!Exists(b_i)) b_i = Sqrt(Abs(a_i^2-Abs(a_m^2-b_m^2))); EndIf

If (!Exists(Nint)) Nint=101; EndIf // No. of nodes on sign-changing interface
If (!Exists(GeomProgint)) GeomProgint=1.1; EndIf // Geometric progression of element size along interface
If (!Exists(Nmu)) Nmu = 4; EndIf // No. of nodes on interfaces [Gamma_I - Gamma] and [Gamma - Gamma_O]
If (!Exists(Nz)) Nz = 10; EndIf // total number of element in z direction
    // Adimensional characteristic length
If (!Exists(CharLengthMin_adim)) CharLengthMin_adim=1; EndIf
If (!Exists(CharLengthMax_adim)) CharLengthMax_adim=2; EndIf

    /**** Quad mesh parameters ****/
If (!Exists(GenerateQuadMesh)) GenerateQuadMesh=1; EndIf
// ---- Preliminary computations

    // Corner parameter (C1 junction)
    // top junction point (ym>0)
x_m = Sqrt(Tan(phi/2)^2+(b_m/a_m)^2);
x_m = -a_m * Tan(phi/2) / x_m;
y_m = b_m*Sqrt(1-(x_m/a_m)^2);
    // abscissa of top point of corner
x_c = x_m - y_m/Tan(phi/2);
    // radius of corner circle
R = Sqrt((x_c-x_m)^2+(y_m)^2); 
R_PML = R_PML*R;
R_TR = R_TR*R;

Printf("Info    : Ellipse gamma-m: (a,b)=(%g,%g) aspect ratio c=%g",a_m,b_m,Sqrt(Abs(a_m^2-b_m^2)));
Printf("Info    : Ellipse gamma-d: (a,b)=(%g,%g) aspect ratio c=%g",a_d,b_d,Sqrt(Abs(a_d^2-b_d^2)));
Printf("Info    : Ellipse gamma-o: (a,b)=(%g,%g) aspect ratio c=%g",a_o,b_o,Sqrt(Abs(a_o^2-b_o^2)));
Printf("Info    : Ellipse gamma-i: (a,b)=(%g,%g) aspect ratio c=%g",a_i,b_i,Sqrt(Abs(a_i^2-b_i^2)));

a = a_i; b = b_i;
x_i = (2*x_c + Sqrt((2*x_c)^2-4*(1-b^2/a^2)*(x_c^2+b^2-R^2)))/(2*(1-b^2/a^2));
y_i = b*Sqrt(1-x_i^2/a^2);
Printf("GammaI top junction point: (%g,%g)",x_i,y_i);
phi_i = 2*Atan(y_i/(x_i-x_c));
Printf("Inner corner angle: %g deg",phi_i*180/Pi);

a = a_o; b = b_o;
x_o = (2*x_c + Sqrt((2*x_c)^2-4*(1-b^2/a^2)*(x_c^2+b^2-R^2)))/(2*(1-b^2/a^2));
y_o = b*Sqrt(1-x_o^2/a^2);
Printf("GammaO top junction point: (%g,%g)",x_o,y_o);
phi_o = 2*Atan(y_o/(x_o-x_c));
Printf("Outer corner angle: %g deg",phi_o*180/Pi);

// Test that R_PML and R_TR are smaller than R
If ((R_PML>R) || (R_TR>R))
    Error("Either R_PML or R_TR is greater than R");
EndIf
// Test for collision with Gamma-D
If (Abs(x_c-R)>a_d)
    Error("Corner circle collides with Gamma-D: reduce R");
EndIf
// check outer ellipse
If ((a_o>a_d) || (a_o<a_m) || (b_o>b_d) || (b_o<b_m))
    Error("Outer ellipse invalid.");
EndIf    
// check inner ellipse
If ((a_i>a_m) || (b_i>b_m))
    Error("Inner ellipse invalid.");
EndIf

    // ---- Geometry
    // Syntax: Ellipse(StartPoint,CenterPoint,a point anywhere on the major axis,EndPoint)
pOrigin = newp; Point(pOrigin) = {0,0,0};

    // Ellipse gamma-d (D)
p_D_Right = newp; Point(p_D_Right) = {a_d,0,0};
p_D_Top = newp; Point(p_D_Top) = {0,b_d,0};
p_D_Left = newp; Point(p_D_Left) = {-a_d,0,0};
p_D_Bot = newp; Point(p_D_Bot) = {0,-b_d,0};
l_D_TR = newl; Ellipse(l_D_TR) = {p_D_Right,pOrigin,p_D_Right,p_D_Top};
l_D_TL = newl; Ellipse(l_D_TL) = {p_D_Top,pOrigin,p_D_Right,p_D_Left};
l_D_BL = newl; Ellipse(l_D_BL) = {p_D_Left,pOrigin,p_D_Right,p_D_Bot};
l_D_BR = newl; Ellipse(l_D_BR) = {p_D_Bot,pOrigin,p_D_Right,p_D_Right};

    // Ellipse gamma-o (O)
p_O_Right = newp; Point(p_O_Right) = {a_o,0,0};
p_O_Top = newp; Point(p_O_Top) = {0,b_o,0};
p_O_LeftTop = newp; Point(p_O_LeftTop) = {x_o,y_o,0};
p_O_LeftBot = newp; Point(p_O_LeftBot) = {x_o,-y_o,0};
p_O_Bottom = newp; Point(p_O_Bottom) = {0,-b_o,0};
l_O_TR = newl; Ellipse(l_O_TR) = {p_O_Right,pOrigin,p_O_Right,p_O_Top};
l_O_TL = newl; Ellipse(l_O_TL) = {p_O_Top,pOrigin,p_O_Right,p_O_LeftTop};
l_O_BL = newl; Ellipse(l_O_BL) = {p_O_LeftBot,pOrigin,p_O_Right,p_O_Bottom};
l_O_BR = newl; Ellipse(l_O_BR) = {p_O_Bottom,pOrigin,p_O_Right,p_O_Right};

    // Ellipse gamma-m (M)
p_M_Right = newp; Point(p_M_Right) = {a_m,0,0};
p_M_Top = newp; Point(p_M_Top) = {0,b_m,0};
p_M_LeftTop = newp; Point(p_M_LeftTop) = {x_m,y_m,0};
p_M_LeftBot = newp; Point(p_M_LeftBot) = {x_m,-y_m,0};
p_M_Bottom = newp; Point(p_M_Bottom) = {0,-b_m,0};
l_M_TR = newl; Ellipse(l_M_TR) = {p_M_Right,pOrigin,p_M_Right,p_M_Top};
l_M_TL = newl; Ellipse(l_M_TL) = {p_M_Top,pOrigin,p_M_Right,p_M_LeftTop};
l_M_BL = newl; Ellipse(l_M_BL) = {p_M_LeftBot,pOrigin,p_M_Right,p_M_Bottom};
l_M_BR = newl; Ellipse(l_M_BR) = {p_M_Bottom,pOrigin,p_M_Right,p_M_Right};

    // Ellipse gamma-i (I)
        /// GammaI interface
p_I_Right = newp; Point(p_I_Right) = {a_i,0,0};
p_I_Top = newp; Point(p_I_Top) = {0,b_i,0};
p_I_LeftTop = newp; Point(p_I_LeftTop) = {x_i,y_i,0};
p_I_LeftBot = newp; Point(p_I_LeftBot) = {x_i,-y_i,0};
p_I_Bottom = newp; Point(p_I_Bottom) = {0,-b_i,0};
l_I_TR = newl; Ellipse(l_I_TR) = {p_I_Right,pOrigin,p_I_Right,p_I_Top};
l_I_TL = newl; Ellipse(l_I_TL) = {p_I_Top,pOrigin,p_I_Right,p_I_LeftTop};
l_I_BL = newl; Ellipse(l_I_BL) = {p_I_LeftBot,pOrigin,p_I_Right,p_I_Bottom};
l_I_BR = newl; Ellipse(l_I_BR) = {p_I_Bottom,pOrigin,p_I_Right,p_I_Right};
    // Corner radius r=R
pC_center = newp; Point(pC_center) = {x_c,0,0};
pC_Left = newp; Point(pC_Left) = {x_c-R,0,0};
        // Outer line (counterclockwise)
lC_Outer_1 = newl; Circle(lC_Outer_1) = {p_I_LeftBot,pC_center,p_I_LeftTop}; 
lC_Outer_2 = newl; Circle(lC_Outer_2) = {p_I_LeftTop,pC_center,p_M_LeftTop}; 
lC_Outer_3 = newl; Circle(lC_Outer_3) = {p_M_LeftTop,pC_center,p_O_LeftTop}; 
lC_Outer_4 = newl; Circle(lC_Outer_4) = {p_O_LeftTop,pC_center,pC_Left}; 
lC_Outer_5 = newl; Circle(lC_Outer_5) = {pC_Left,pC_center,p_O_LeftBot}; 
lC_Outer_6 = newl; Circle(lC_Outer_6) = {p_O_LeftBot,pC_center,p_M_LeftBot}; 
lC_Outer_7 = newl; Circle(lC_Outer_7) = {p_M_LeftBot,pC_center,p_I_LeftBot}; 


    // Domain omega-o (outside corner disc)
l_OmegaO_Right = newl; Line(l_OmegaO_Right)={p_M_Right,p_O_Right};
l_OmegaO_Top = newl; Line(l_OmegaO_Top)={p_M_Top,p_O_Top};
l_OmegaO_Bottom = newl; Line(l_OmegaO_Bottom)={p_M_Bottom,p_O_Bottom};

l_OmegaO_TR = newl; Line Loop(l_OmegaO_TR) = {l_OmegaO_Right,l_O_TR,-l_OmegaO_Top,-l_M_TR};
s_OmegaO_TR = news; Plane Surface(s_OmegaO_TR)={l_OmegaO_TR};

l_OmegaO_BR = newl; Line Loop(l_OmegaO_BR) = {l_OmegaO_Bottom,l_O_BR,-l_OmegaO_Right,-l_M_BR};
s_OmegaO_BR = news; Plane Surface(s_OmegaO_BR)={l_OmegaO_BR};

l_OmegaO_TL = newl; Line Loop(l_OmegaO_TL) = {l_OmegaO_Top,l_O_TL,-lC_Outer_3,-l_M_TL};
s_OmegaO_TL = news; Plane Surface(s_OmegaO_TL)={l_OmegaO_TL};

l_OmegaO_BL = newl; Line Loop(l_OmegaO_BL) = {-lC_Outer_6,l_O_BL,-l_OmegaO_Bottom,-l_M_BL};
s_OmegaO_BL = news; Plane Surface(s_OmegaO_BL)={l_OmegaO_BL};

    // Domain omega-i (outside corner disc)
l_OmegaI_Right = newl; Line(l_OmegaI_Right)={p_I_Right,p_M_Right};
l_OmegaI_Top = newl; Line(l_OmegaI_Top)={p_I_Top,p_M_Top};
l_OmegaI_Bottom = newl; Line(l_OmegaI_Bottom)={p_I_Bottom,p_M_Bottom};

l_OmegaI_TR = newl; Line Loop(l_OmegaI_TR) = {l_OmegaI_Right,l_M_TR,-l_OmegaI_Top,-l_I_TR};
s_OmegaI_TR = news; Plane Surface(s_OmegaI_TR)={l_OmegaI_TR};

l_OmegaI_TL = newl; Line Loop(l_OmegaI_TL) = {l_OmegaI_Top,l_M_TL,-lC_Outer_2,-l_I_TL};
s_OmegaI_TL = news; Plane Surface(s_OmegaI_TL)={l_OmegaI_TL};

l_OmegaI_BL = newl; Line Loop(l_OmegaI_BL) = {-lC_Outer_7,l_M_BL,-l_OmegaI_Bottom,-l_I_BL};
s_OmegaI_BL = news; Plane Surface(s_OmegaI_BL)={l_OmegaI_BL};

l_OmegaI_BR = newl; Line Loop(l_OmegaI_BR) = {l_OmegaI_Bottom,l_M_BR,-l_OmegaI_Right,-l_I_BR};
s_OmegaI_BR = news; Plane Surface(s_OmegaI_BR)={l_OmegaI_BR};

    // ---- Mesh (structured region)
        // Orthoradial distribution on circle {r=R}
Transfinite Line{lC_Outer_2,lC_Outer_3,lC_Outer_6,lC_Outer_7} = Nmu;
// Transfinite Line{lC_Outer_1} = 2*Nmu*Floor(phi/(phi_o/2-phi_i/2));
// Transfinite Line{lC_Outer_4,lC_Outer_5} = 2*Nmu*Floor((Pi-phi_o/2)/(phi_o/2-phi_i/2));

//         Structured mesh on omega-o (outside corner disc)
Transfinite Line{l_OmegaO_Right,l_OmegaO_Top,l_OmegaO_Bottom} = Nmu;
Transfinite Line{l_O_TR,l_M_TR} = Floor(Nint/4) Using Progression GeomProgint;
Transfinite Line{-l_O_TL,-l_M_TL} = Floor(Nint/4) Using Progression GeomProgint;
Transfinite Line{l_O_BL,l_M_BL} = Floor(Nint/4) Using Progression GeomProgint;
Transfinite Line{-l_O_BR,-l_M_BR} = Floor(Nint/4) Using Progression GeomProgint;

Transfinite Surface(s_OmegaO_TR) = {p_M_Right,p_O_Right,p_M_Top,p_O_Top} Alternate;
Transfinite Surface(s_OmegaO_TL) = {p_M_LeftTop,p_O_LeftTop,p_M_Top,p_O_Top} Alternate;
Transfinite Surface(s_OmegaO_BL) = {p_M_LeftBot,p_O_LeftBot,p_M_Bottom,p_O_Bottom} Alternate;
Transfinite Surface(s_OmegaO_BR) = {p_M_Right,p_O_Right,p_M_Bottom,p_O_Bottom} Alternate;

//         Structured mesh on omega-i (outside corner disc)
Transfinite Line{l_OmegaI_Right,l_OmegaI_Top,l_OmegaI_Bottom} = Nmu;
Transfinite Line{l_I_TR,-l_I_TL,l_I_BL,-l_I_BR} = Floor(Nint/4) Using Progression GeomProgint;

Transfinite Surface(s_OmegaI_TR) = {p_I_Right,p_M_Right,p_M_Top,p_I_Top} Alternate;
Transfinite Surface(s_OmegaI_TL) = {p_I_LeftTop,p_M_LeftTop,p_M_Top,p_I_Top} Alternate;
Transfinite Surface(s_OmegaI_BL) = {p_I_LeftBot,p_M_LeftBot,p_M_Bottom,p_I_Bottom} Alternate;
Transfinite Surface(s_OmegaI_BR) = {p_I_Right,p_M_Right,p_M_Bottom,p_I_Bottom} Alternate;

    // ---- Mesh (unstructured region)
    /* Structured interface controlled by Nint and Ny.
     * Unstructured part is coarsened away from the interface.
     */

    // Diagnostic messages
perimeter_GammaI = 2*Pi*Sqrt(a_m*b_m);    // approximate value (lower bound)
Printf("Gamma_I perimeter: %g",perimeter_GammaI);
lc = perimeter_GammaI/(Nint); // element size in structured region (lower bound)
Printf("Characteristic length: %g",lc);

// Mesh.CharacteristicLengthMin = lc*CharLengthMin_adim;
// Mesh.CharacteristicLengthMax = lc*CharLengthMax_adim;

// When the element size is fully specified by a background mesh field, set:
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;

    // Unstructured part controlled by Distance and Threshold fields
Printf("Distance+Threshold params...");
Field[1]=Distance;
Field[1].EdgesList = {l_O_BL,l_O_BR,l_O_TR,l_O_TL,l_I_BL,l_I_BR,l_I_TR,l_I_TL,lC_Outer_1,lC_Outer_4,lC_Outer_5};
Field[1].NNodesByEdge = 100;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = lc*CharLengthMin_adim;
Field[2].LcMax = lc*CharLengthMax_adim;
Field[2].DistMin = 0;
Field[2].DistMax = 0.5*(a_d-a_o);
Field[2].Sigmoid=1;
Background Field = 2;

    // ---- Corner geometry
name = "corner";
    // Abscissa
x_input[] = {Log(R_TR), Log(R_PML), Log(R)};
    // Angle
y_input[] = {-Pi, -phi_o/2, -phi/2,-phi_i/2,phi_i/2,phi/2,phi_o/2,Pi};

For i In {0:(#x_input[]-1)}
    x_input[i] = x_input[i] + x_ofst;
EndFor

For i In {0:(#y_input[]-1)}
    y_input[i] = y_input[i] + y_ofst;
EndFor

    // Mesh in z-direction
For i In {0:(#x_input[]-2)}
    N_mesh_x[i] = Ceil(Nz*(x_input[i+1]-x_input[i])/(x_input[#x_input[]-1]-x_input[0]))+1;
EndFor

    // Definition of physical entities
    // Set to '0' to disable default physical entities
Bool_define_line_physical_entities = 1;
    // disable default surface physical entities
Bool_define_surface_physical_entities = 0;

Include "Macro-Rectangle-Array.geo";

     // Periodicity disk boundary -> Right boundary
l_src[] = {lC_Outer_5, lC_Outer_6, lC_Outer_7, lC_Outer_1, lC_Outer_2, lC_Outer_3, lC_Outer_4};
Call ConnectSrcToRightBnd;

    // ----

    // Domains omega-d and omega-m (outside corner disk)
    // Recombination algorithm requires non-transfinite regions
    // to be meshed last.
l_D = newl; Line Loop(l_D) = {l_D_BL,l_D_BR,l_D_TR,l_D_TL};
l_O = newl; Line Loop(l_O) = {l_O_BL,l_O_BR,l_O_TR,l_O_TL,lC_Outer_4,lC_Outer_5};
l_I = newl; Line Loop(l_I) = {l_I_BL,l_I_BR,l_I_TR,l_I_TL,-lC_Outer_1};

s_OmegaD = news; Plane Surface(s_OmegaD) = {l_D,l_O};
s_OmegaM = news; Plane Surface(s_OmegaM) = {l_I};

    // ---- Physical entities
Physical Line("gamma-d") = {l_D_TR,l_D_TL,l_D_BL,l_D_BR};
Physical Line("gamma-m") = {l_M_TR,l_M_TL,l_M_BL,l_M_BR};
Physical Surface("omega-m") = {s_OmegaM,s_OmegaI_TR,s_OmegaI_TL,s_OmegaI_BL,s_OmegaI_BR};
Physical Surface("omega-d") = {s_OmegaD,s_OmegaO_TR,s_OmegaO_TL,s_OmegaO_BL,s_OmegaO_BR};

    // For convenience, define new physical entities
Physical Surface("corner-omega-m") = {s~{7},s~{10},s~{13}};
Physical Surface("corner-omega-d") = {s~{1},s~{4},s~{16},s~{19}};
Physical Surface("corner-omega-d-pml") = {s~{0},s~{3},s~{15},s~{18}};
Physical Surface("corner-omega-m-pml") = {s~{6},s~{9},s~{12}};

    // -- Meshing

//Mesh.RecombinationAlgorithm = 3;
//Mesh.Algorithm = 6;
//Mesh.ElementOrder = 2;
If (GenerateQuadMesh>0)
    Include 'Macro-Generate-Quad.geo';
Else
    Mesh 2;
EndIf

/*
    // Recombination algorithm are very sensitive to mesh size.
    // If recombination fails or is only partial (i.e. there are
    // triangles remaining) then tweak mesh size at interface and/or
    // characteristic mesh size.
If (GenerateQuadMesh)
    Mesh.RecombineAll = 1;
    // Mesh recombination algorithm (0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad)
    Mesh.RecombinationAlgorithm = 3;
    // 2D mesh algorithm (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms)
    Mesh.Algorithm = 9;
EndIf
*/
