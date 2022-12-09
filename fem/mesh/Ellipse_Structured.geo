/*********************************************************************
 * Unstructured mesh of an ellipse, structured mesh at interface.
 *
 *  - Mesh structured at sign-changing interface, defined using two
 * intermediary ellipses (denoted 'i' for inner and 'o' for outer).
 *  - Recombination algorithm can be used to produce quad mesh.
 * /

    /*****************************************************************
     * Input parameters
     * **************************************************************/ 
        /**** Geometrical parameters ****/

            /* Ellipse gamma-m
            *    a_m, b_m */
If (!Exists(a_m)) a_m=2.5; EndIf
If (!Exists(b_m)) b_m=1; EndIf
c = Sqrt(Abs(a_m^2-b_m^2));

            /* Ellipse gamma-d
            *    a_d, b_d        semi-axes */
If (!Exists(a_d)) a_d=3; EndIf
If (!Exists(b_d)) b_d=Sqrt(a_d^2-c^2); EndIf

        // Ellipse gamma-o
If (!Exists(a_o)) a_o=1.02*a_m; EndIf
If (!Exists(b_o)) b_o=Sqrt(a_o^2-c^2); EndIf

                // Ellipse gamma-i
If (!Exists(a_i)) a_i=0.98*a_m; EndIf
If (!Exists(b_i)) b_i=Sqrt(a_i^2-c^2); EndIf

    /**** Mesh parameters (Structured interface) ****/
If (!Exists(N_mu)) N_mu=1; EndIf
If (!Exists(Nint)) Nint=50; EndIf
If (!Exists(GeomProgint)) GeomProgint=1.01; EndIf
If (!Exists(TransfiniteSurface)) TransfiniteSurface=1; EndIf
    // Adimensional characteristic length
If (!Exists(CharLengthMin_adim)) CharLengthMin_adim=1; EndIf
If (!Exists(CharLengthMax_adim)) CharLengthMax_adim=5; EndIf

    /**** Quad mesh parameters ****/
If (!Exists(GenerateQuadMesh)) GenerateQuadMesh=0; EndIf

pOrigin = newp; Point(pOrigin) = {0,0,0};
    // Ellipse(StartPoint,CenterPoint,a point anywhere on the major axis,EndPoint)
    // Ellipse gamma-d
c = Sqrt(Abs(a_d^2-b_d^2));
Printf("Info    : Ellipse gamma-d: (a,b)=(%g,%g) aspect ratio c=%g",a_d,b_d,c);
p_GammaD_Right = newp; Point(p_GammaD_Right) = {a_d,0,0};
p_GammaD_Top = newp; Point(p_GammaD_Top) = {0,b_d,0};
p_GammaD_Left = newp; Point(p_GammaD_Left) = {-a_d,0,0};
p_GammaD_Bottom = newp; Point(p_GammaD_Bottom) = {0,-b_d,0};
l_GammaD_TR = newl; Ellipse(l_GammaD_TR) = {p_GammaD_Right,pOrigin,pOrigin,p_GammaD_Top};
l_GammaD_TL = newl; Ellipse(l_GammaD_TL) = {p_GammaD_Top,pOrigin,pOrigin,p_GammaD_Left};
l_GammaD_BL = newl; Ellipse(l_GammaD_BL) = {p_GammaD_Left,pOrigin,pOrigin,p_GammaD_Bottom};
l_GammaD_BR = newl; Ellipse(l_GammaD_BR) = {p_GammaD_Bottom,pOrigin,pOrigin,p_GammaD_Right};

    // Ellipse gamma-o
Printf("Info    : Ellipse gamma-o: (a,b)=(%g,%g) aspect ratio c=%g",a_o,b_o,c);
p_GammaO_Right = newp; Point(p_GammaO_Right) = {a_o,0,0};
p_GammaO_Top = newp; Point(p_GammaO_Top) = {0,b_o,0};
p_GammaO_Left = newp; Point(p_GammaO_Left) = {-a_o,0,0};
p_GammaO_Bottom = newp; Point(p_GammaO_Bottom) = {0,-b_o,0};
l_GammaO_TR = newl; Ellipse(l_GammaO_TR) = {p_GammaO_Right,pOrigin,pOrigin,p_GammaO_Top};
l_GammaO_TL = newl; Ellipse(l_GammaO_TL) = {p_GammaO_Top,pOrigin,pOrigin,p_GammaO_Left};
l_GammaO_BL = newl; Ellipse(l_GammaO_BL) = {p_GammaO_Left,pOrigin,pOrigin,p_GammaO_Bottom};
l_GammaO_BR = newl; Ellipse(l_GammaO_BR) = {p_GammaO_Bottom,pOrigin,pOrigin,p_GammaO_Right};

    // Ellipse gamma-m
Printf("Info    : Ellipse gamma-m: (a,b)=(%g,%g) aspect ratio c=%g",a_m,b_m,c);
p_GammaM_Right = newp; Point(p_GammaM_Right) = {a_m,0,0};
p_GammaM_Top = newp; Point(p_GammaM_Top) = {0,b_m,0};
p_GammaM_Left = newp; Point(p_GammaM_Left) = {-a_m,0,0};
p_GammaM_Bottom = newp; Point(p_GammaM_Bottom) = {0,-b_m,0};
l_GammaM_TR = newl; Ellipse(l_GammaM_TR) = {p_GammaM_Right,pOrigin,pOrigin,p_GammaM_Top};
l_GammaM_TL = newl; Ellipse(l_GammaM_TL) = {p_GammaM_Top,pOrigin,pOrigin,p_GammaM_Left};
l_GammaM_BL = newl; Ellipse(l_GammaM_BL) = {p_GammaM_Left,pOrigin,pOrigin,p_GammaM_Bottom};
l_GammaM_BR = newl; Ellipse(l_GammaM_BR) = {p_GammaM_Bottom,pOrigin,pOrigin,p_GammaM_Right};

    // Ellipse gamma-i
Printf("Info    : Ellipse gamma-i: (a,b)=(%g,%g) aspect ratio c=%g",a_i,b_i,c);
p_GammaI_Right = newp; Point(p_GammaI_Right) = {a_i,0,0};
p_GammaI_Top = newp; Point(p_GammaI_Top) = {0,b_i,0};
p_GammaI_Left = newp; Point(p_GammaI_Left) = {-a_i,0,0};
p_GammaI_Bottom = newp; Point(p_GammaI_Bottom) = {0,-b_i,0};
l_GammaI_TR = newl; Ellipse(l_GammaI_TR) = {p_GammaI_Right,pOrigin,pOrigin,p_GammaI_Top};
l_GammaI_TL = newl; Ellipse(l_GammaI_TL) = {p_GammaI_Top,pOrigin,pOrigin,p_GammaI_Left};
l_GammaI_BL = newl; Ellipse(l_GammaI_BL) = {p_GammaI_Left,pOrigin,pOrigin,p_GammaI_Bottom};
l_GammaI_BR = newl; Ellipse(l_GammaI_BR) = {p_GammaI_Bottom,pOrigin,pOrigin,p_GammaI_Right};


    // Domain omega-o
l_OmegaO_Right = newl; Line(l_OmegaO_Right)={p_GammaM_Right,p_GammaO_Right};
l_OmegaO_Top = newl; Line(l_OmegaO_Top)={p_GammaM_Top,p_GammaO_Top};
l_OmegaO_Left = newl; Line(l_OmegaO_Left)={p_GammaM_Left,p_GammaO_Left};
l_OmegaO_Bottom = newl; Line(l_OmegaO_Bottom)={p_GammaM_Bottom,p_GammaO_Bottom};

l_OmegaO_TR = newl; Line Loop(l_OmegaO_TR) = {l_OmegaO_Right,l_GammaO_TR,-l_OmegaO_Top,-l_GammaM_TR};
s_OmegaO_TR = news; Plane Surface(s_OmegaO_TR)={l_OmegaO_TR};

l_OmegaO_TL = newl; Line Loop(l_OmegaO_TL) = {l_OmegaO_Top,l_GammaO_TL,-l_OmegaO_Left,-l_GammaM_TL};
s_OmegaO_TL = news; Plane Surface(s_OmegaO_TL)={l_OmegaO_TL};

l_OmegaO_BL = newl; Line Loop(l_OmegaO_BL) = {l_OmegaO_Left,l_GammaO_BL,-l_OmegaO_Bottom,-l_GammaM_BL};
s_OmegaO_BL = news; Plane Surface(s_OmegaO_BL)={l_OmegaO_BL};

l_OmegaO_BR = newl; Line Loop(l_OmegaO_BR) = {l_OmegaO_Bottom,l_GammaO_BR,-l_OmegaO_Right,-l_GammaM_BR};
s_OmegaO_BR = news; Plane Surface(s_OmegaO_BR)={l_OmegaO_BR};


    // Domain omega-i
l_OmegaI_Right = newl; Line(l_OmegaI_Right)={p_GammaI_Right,p_GammaM_Right};
l_OmegaI_Top = newl; Line(l_OmegaI_Top)={p_GammaI_Top,p_GammaM_Top};
l_OmegaI_Left = newl; Line(l_OmegaI_Left)={p_GammaI_Left,p_GammaM_Left};
l_OmegaI_Bottom = newl; Line(l_OmegaI_Bottom)={p_GammaI_Bottom,p_GammaM_Bottom};

l_OmegaI_TR = newl; Line Loop(l_OmegaI_TR) = {l_OmegaI_Right,l_GammaM_TR,-l_OmegaI_Top,-l_GammaI_TR};
s_OmegaI_TR = news; Plane Surface(s_OmegaI_TR)={l_OmegaI_TR};

l_OmegaI_TL = newl; Line Loop(l_OmegaI_TL) = {l_OmegaI_Top,l_GammaM_TL,-l_OmegaI_Left,-l_GammaI_TL};
s_OmegaI_TL = news; Plane Surface(s_OmegaI_TL)={l_OmegaI_TL};

l_OmegaI_BL = newl; Line Loop(l_OmegaI_BL) = {l_OmegaI_Left,l_GammaM_BL,-l_OmegaI_Bottom,-l_GammaI_BL};
s_OmegaI_BL = news; Plane Surface(s_OmegaI_BL)={l_OmegaI_BL};

l_OmegaI_BR = newl; Line Loop(l_OmegaI_BR) = {l_OmegaI_Bottom,l_GammaM_BR,-l_OmegaI_Right,-l_GammaI_BR};
s_OmegaI_BR = news; Plane Surface(s_OmegaI_BR)={l_OmegaI_BR};

        // Structured mesh on omega-o
Transfinite Line{l_OmegaO_Right,l_OmegaO_Top,l_OmegaO_Left,l_OmegaO_Bottom} = N_mu+1;
Transfinite Line{l_GammaO_TR,-l_GammaO_TL,l_GammaO_BL,-l_GammaO_BR} = Nint+1 Using Progression GeomProgint;
Transfinite Line{l_GammaM_TR,-l_GammaM_TL,l_GammaM_BL,-l_GammaM_BR} = Nint+1 Using Progression GeomProgint;
If (TransfiniteSurface)
    Transfinite Surface(s_OmegaO_TR) = {p_GammaM_Right,p_GammaO_Right,p_GammaM_Top,p_GammaO_Top} AlternateLeft;
    Transfinite Surface(s_OmegaO_TL) = {p_GammaM_Left,p_GammaO_Left,p_GammaM_Top,p_GammaO_Top} Alternate;
    Transfinite Surface(s_OmegaO_BL) = {p_GammaM_Left,p_GammaO_Left,p_GammaM_Bottom,p_GammaO_Bottom} Alternate;
    Transfinite Surface(s_OmegaO_BR) = {p_GammaM_Right,p_GammaO_Right,p_GammaM_Bottom,p_GammaO_Bottom} Alternate;
EndIf
// 
// Recombine Surface {s_OmegaO_TR,s_OmegaO_TL,s_OmegaO_BL,s_OmegaO_BR};

        // Structured mesh on omega-I
Transfinite Line{l_OmegaI_Right,l_OmegaI_Top,l_OmegaI_Left,l_OmegaI_Bottom} = N_mu+1;
Transfinite Line{l_GammaI_TR,-l_GammaI_TL,l_GammaI_BL,-l_GammaI_BR} = Nint+1 Using Progression GeomProgint;
If (TransfiniteSurface)
    Transfinite Surface(s_OmegaI_TR) = {p_GammaI_Right,p_GammaM_Right,p_GammaM_Top,p_GammaI_Top} Alternate;
    Transfinite Surface(s_OmegaI_TL) = {p_GammaI_Left,p_GammaM_Left,p_GammaM_Top,p_GammaI_Top} Alternate;
    Transfinite Surface(s_OmegaI_BL) = {p_GammaI_Left,p_GammaM_Left,p_GammaM_Bottom,p_GammaI_Bottom} Alternate;
    Transfinite Surface(s_OmegaI_BR) = {p_GammaI_Right,p_GammaM_Right,p_GammaM_Bottom,p_GammaI_Bottom} Alternate;
EndIf
// 
// Recombine Surface {s_OmegaI_TR,s_OmegaI_TL,s_OmegaI_BL,s_OmegaI_BR};

    // Domain omega-d
    // Recombination algorithm requires non-transfinite surfaces
    // to be meshed last.
l_OmegaD_Right = newl; Line(l_OmegaD_Right)={p_GammaO_Right,p_GammaD_Right};
l_OmegaD_Top = newl; Line(l_OmegaD_Top)={p_GammaO_Top,p_GammaD_Top};
l_OmegaD_Left = newl; Line(l_OmegaD_Left)={p_GammaO_Left,p_GammaD_Left};
l_OmegaD_Bottom = newl; Line(l_OmegaD_Bottom)={p_GammaO_Bottom,p_GammaD_Bottom};

l_OmegaD_TR = newl; Line Loop(l_OmegaD_TR) = {l_OmegaD_Right,l_GammaD_TR,-l_OmegaD_Top,-l_GammaO_TR};
s_OmegaD_TR = news; Plane Surface(s_OmegaD_TR)={l_OmegaD_TR};

l_OmegaD_TL = newl; Line Loop(l_OmegaD_TL) = {l_OmegaD_Top,l_GammaD_TL,-l_OmegaD_Left,-l_GammaO_TL};
s_OmegaD_TL = news; Plane Surface(s_OmegaD_TL)={l_OmegaD_TL};

l_OmegaD_BL = newl; Line Loop(l_OmegaD_BL) = {l_OmegaD_Left,l_GammaD_BL,-l_OmegaD_Bottom,-l_GammaO_BL};
s_OmegaD_BL = news; Plane Surface(s_OmegaD_BL)={l_OmegaD_BL};

l_OmegaD_BR = newl; Line Loop(l_OmegaD_BR) = {l_OmegaD_Bottom,l_GammaD_BR,-l_OmegaD_Right,-l_GammaO_BR};
s_OmegaD_BR = news; Plane Surface(s_OmegaD_BR)={l_OmegaD_BR};

    // Domain omega-m
l_OmegaM_Right = newl; Line(l_OmegaM_Right)={pOrigin,p_GammaI_Right};
l_OmegaM_Top = newl; Line(l_OmegaM_Top)={pOrigin,p_GammaI_Top};
l_OmegaM_Left = newl; Line(l_OmegaM_Left)={pOrigin,p_GammaI_Left};
l_OmegaM_Bottom = newl; Line(l_OmegaM_Bottom)={pOrigin,p_GammaI_Bottom};

l_OmegaM_TR = newl; Line Loop(l_OmegaM_TR) = {l_OmegaM_Right,l_GammaI_TR,-l_OmegaM_Top};
s_OmegaM_TR = news; Plane Surface(s_OmegaM_TR)={l_OmegaM_TR};

l_OmegaM_TL = newl; Line Loop(l_OmegaM_TL) = {l_OmegaM_Left,-l_GammaI_TL,-l_OmegaM_Top};
s_OmegaM_TL = news; Plane Surface(s_OmegaM_TL)={l_OmegaM_TL};

l_OmegaM_BL = newl; Line Loop(l_OmegaM_BL) = {l_OmegaM_Left,l_GammaI_BL,-l_OmegaM_Bottom};
s_OmegaM_BL = news; Plane Surface(s_OmegaM_BL)={l_OmegaM_BL};

l_OmegaM_BR = newl; Line Loop(l_OmegaM_BR) = {l_OmegaM_Right,-l_GammaI_BR,-l_OmegaM_Bottom};
s_OmegaM_BR = news; Plane Surface(s_OmegaM_BR)={l_OmegaM_BR};

    // Transfinite distribution on x-axis for smoother mesh
//Transfinite Line{-l_OmegaM_Right,-l_OmegaM_Left} = 15 Using Progression 1.20;
    // Transfinite distribution on y-axis for smoother mesh
//Transfinite Line{-l_OmegaM_Top,-l_OmegaM_Bottom} = 15 Using Progression 1.20;

    /************************************************************************
     * Definition of physical elements
     * !! For XLiFE, no space in domain names.
     * **********************************************************************/

Physical Line("gamma-d") = {l_GammaD_TR,l_GammaD_TL,l_GammaD_BL,l_GammaD_BR};
Physical Line("gamma-m") = {l_GammaM_TR,l_GammaM_TL,l_GammaM_BL,l_GammaM_BR};
Physical Surface("omega-m") = {s_OmegaM_TR,s_OmegaI_TR,s_OmegaM_TL,s_OmegaI_TL,s_OmegaM_BL,s_OmegaI_BL,s_OmegaM_BR,s_OmegaI_BR};
Physical Surface("omega-d") = {s_OmegaO_TR,s_OmegaD_TR,s_OmegaO_TL,s_OmegaD_TL,s_OmegaO_BL,s_OmegaD_BL,s_OmegaO_BR,s_OmegaD_BR};

    /************************************************************************
     * Mesh size definition
     * Structured interface controlled by Nint and N_mu.
     * Unstructured part is coarsened away from the interface.
     * **********************************************************************/

    // Diagnostic messages
perimeter_GammaI = 2*Pi*Sqrt(a_m*b_m);    // approximate value (lower bound)
Printf("Gamma_I perimeter: %g",perimeter_GammaI); 
lc = perimeter_GammaI/(4*Nint); // element size in structured region (lower bound)
Printf("Characteristic length: %g",lc);

		/* Uniform unstructured mesh */
// Mesh.CharacteristicLengthMin = CharLengthMin_adim*lc;
// Mesh.CharacteristicLengthMax = CharLengthMax_adim*lc;


// When the element size is fully specified by a background mesh field, set:
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthExtendFromBoundary = 0;


//     Unstructured part controlled by Distance and Threshold fields
Printf("Distance+Threshold params...");
Field[1]=Distance;
Field[1].EdgesList = {l_GammaO_TR,l_GammaI_TR,l_GammaO_TL,l_GammaI_TL,l_GammaO_BR,l_GammaI_BR,l_GammaO_BL,l_GammaI_BL};
Field[1].NNodesByEdge = 5*Nint;

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = CharLengthMin_adim*lc;
Field[2].LcMax = CharLengthMax_adim*lc;

If (b_o-b_m>b_m-b_i)
    Printf("dMin = b_o-b_m = %g",b_o-b_m);
    Field[2].DistMin = b_o-b_m;
Else
    Printf("dMin = b_m-b_i = %g",b_o-b_m);
    Field[2].DistMin = b_m-b_i;
EndIf

Field[2].DistMin = 0;
Field[2].DistMax = 0.5*(a_o-a_m);
Field[2].Sigmoid=1;
Background Field = 2;

    // -- Meshing
//Mesh.RecombinationAlgorithm = 3;
//Mesh.Algorithm = 9;
//Mesh.ElementOrder = 2;

If (GenerateQuadMesh>0)
    Include 'Macro-Generate-Quad.geo';
Else
    Mesh 2;
EndIf
