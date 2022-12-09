/*********************************************************************
 * Ellipse perturbed by one corner on the major axis.
 *
 *  - This gmsh file has been partially assembled using the GUI.
 *  - The Euler region is on the disk of radius R.
 *  - Unstructured mesh, structured at the sign-changing interface.
 *  - Corner at (x_c,0), symmetric w.r.t. x-axis
 *  - Corner disk discretized in Euler coordinates (z,theta)
 *  - Recombination algorithm can be used to produce quad mesh.
 */

 /* Possible improvements:
 *   - Adjust element distribution on sign-changing interface (too dense on
 * top left boundary.)
 *   - Implemenent mesh coarsening using background fields.
 */

// BooleanFragments yields a different line numbering in gmsh 4.9.4 and
// newer, which breaks the geometry.
Warning("This .geo file fails with gmsh 4.9.4 and newer");

SetFactory("OpenCASCADE");

    // ---- Input parameters

    // Ellipse gamma-m semi-axes
If (!Exists(a_m)) a_m=2.5; EndIf
If (!Exists(b_m)) b_m=2; EndIf
c = Abs(a_m^2-b_m^2);
    // Ellipse gamma-d semi-axes
If (!Exists(a_d)) a_d=1.7*a_m; EndIf
If (!Exists(b_d)) b_d=Sqrt(a_d^2-c); EndIf
    // Left corner
If (!Exists(phi1)) phi1=0.75*Pi; EndIf // angle in (0,Pi) (rad)
If (!Exists(R1_TR)) R1_TR = 0.05; EndIf // R_TR/R radius of truncated region
    // Corner region offset
If (!Exists(x_ofst_1)) x_ofst_1 = -a_d-0.5; EndIf
If (!Exists(y_ofst_1)) y_ofst_1 = 0; EndIf

    // Mesh parameters (Structured interface)
    // outer ellipse
If (!Exists(a_o)) a_o=(1-0.1)*a_m+0.1*a_d; EndIf
If (!Exists(b_o)) b_o=Sqrt(Abs(a_o^2-c)); EndIf
    // inner ellipse
If (!Exists(a_i)) a_i=(1-5e-2)*a_m; EndIf
If (!Exists(b_i)) b_i = Sqrt(Abs(a_i^2-c)); EndIf

If (!Exists(Nint)) Nint=101; EndIf // No. of nodes on sign-changing interface
If (!Exists(GeomProgint)) GeomProgint=1.01; EndIf // Geometric progression of element size along interface
If (!Exists(Nmu)) Nmu = 1; EndIf // No. of nodes on interfaces [Gamma_I - Gamma] and [Gamma - Gamma_O]
If (!Exists(Nz)) Nz = 10; EndIf // total number of element in z direction
    // Adimensional characteristic length
If (!Exists(CharLengthMin_adim)) CharLengthMin_adim=1; EndIf
If (!Exists(CharLengthMax_adim)) CharLengthMax_adim=2; EndIf
    // Recombination algorithm are very sensitive to mesh size.
    // If recombination fails or is only partial (i.e. there are
    // triangles remaining) then tweak mesh size at interface and/or
    // characteristic mesh size.
    // See 'Unstructured quadrangular meshes' gmsh tutorial.
If (!Exists(GenerateQuadMesh)) GenerateQuadMesh=1; EndIf

    // ---- Ellipses
p=newp; Point(p) = {0, 0, 0, 0};
a=a_i; b=b_i;
p=newl; Ellipse(p) = {0, 0, 0, a, b, 0, 0.5*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a, b, 0.5*Pi, Pi};
p=newl; Ellipse(p) = {0, 0, 0, a, b, Pi, 1.5*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a, b, 1.5*Pi, 2*Pi};
a=a_m; b=b_m;
p=newl; Ellipse(p) = {0, 0, 0, a, b, 0, 0.5*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a, b, 0.5*Pi, Pi};
p=newl; Ellipse(p) = {0, 0, 0, a, b, Pi, 1.5*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a, b, 1.5*Pi, 2*Pi};
a=a_o; b=b_o;
p=newl; Ellipse(p) = {0, 0, 0, a, b, 0, 0.5*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a, b, 0.5*Pi, Pi};
p=newl; Ellipse(p) = {0, 0, 0, a, b, Pi, 1.5*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a, b, 1.5*Pi, 2*Pi};
a=a_d; b=b_d;
p=newl; Ellipse(p) = {0, 0, 0, a, b, 0, 0.5*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a, b, 0.5*Pi, Pi};
p=newl; Ellipse(p) = {0, 0, 0, a, b, Pi, 1.5*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a, b, 1.5*Pi, 2*Pi};

    // ---- Points associated with corner on major axis (C^1 junction)
    // Junction point on gamma-m
    // (slope tan(phi1/2) and yj>0)
x_j = Sqrt(Tan(phi1/2)^2+(b_m/a_m)^2);
x_j = -a_m * Tan(phi1/2) / x_j;
y_j = b_m*Sqrt(1-(x_j/a_m)^2);
    // corner and radius
x_c = x_j - y_j/Tan(phi1/2);
y_c = 0;
p = newp; Point(p) = {x_c,y_c,0};
R = Sqrt((x_c-x_j)^2+(y_c-y_j)^2);
p = newl; Circle(p) = {x_c, y_c, 0, R, 0, Pi};
p = newl; Circle(p) = {x_c, y_c, 0, R, Pi, 2*Pi};
If (Abs(x_c-R)>a_d)
    Error("Corner circle collides with Gamma-D: increase phi1");
EndIf
R1 = R;
    // Junction points
BooleanFragments{ Curve{17}; Curve{18}; Delete; }{ Curve{10}; Curve{6}; Curve{2}; Curve{11}; Curve{7}; Curve{3}; Delete; }
Recursive Delete {
  Curve{26}; Curve{28}; Curve{30}; Curve{35}; Curve{33}; Curve{31};
}

    // Angle of each junction point

xj[] = Point{36}; // gamma-o
phi_o = 2*Atan((xj[1]-y_c)/(xj[0]-x_c));
Printf("[Left corner] Angle outer: %g deg",phi_o*(180/Pi));
xj[] = Point{37}; // gamma-m
phi_m = 2*Atan((xj[1]-y_c)/(xj[0]-x_c));
Printf("[Left corner] Angle m: %g deg",phi_m*(180/Pi));
xj[] = Point{38}; // gamma-i
phi_i = 2*Atan((xj[1]-y_c)/(xj[0]-x_c));
Printf("[Left corner] Angle inner: %g deg",phi_i*(180/Pi));
phi1_o = phi_o; phi1_m = phi_m; phi1_i = phi_i;


    // ---- Lines to subdivise interface boundary layer

Line(37) = {2, 10};
//+
Line(38) = {10, 18};
//+
Line(39) = {3, 11};
//+
Line(40) = {11, 19};
//+
Line(41) = {8, 16};
//+
Line(42) = {16, 24};


    // ---- Left corner geometry in Euler coordinates
    // (A) Define points
    // Abscissa
x_input[] = {Log(R1*R1_TR), Log(R1)};
    // Angle
y_input[] = {-Pi, -phi1_o/2, -phi1_m/2,-phi1_i/2,0,phi1_i/2,phi1_m/2,phi1_o/2,Pi};

For i In {0:(#x_input[]-1)}
    x_input[i] = x_input[i] + x_ofst_1;
EndFor

For i In {0:(#y_input[]-1)}
    y_input[i] = y_input[i] + y_ofst_1;
EndFor

For i In {0:(#x_input[]-1)}
    For j In {0:(#y_input[]-1)}
        p = newp; Point(p) = {x_input[i],y_input[j],0};
    EndFor
EndFor

    // (B) Define lines
Line(43) = {69, 60};
//+
Line(44) = {68, 59};
//+
Line(45) = {67, 58};
//+
Line(46) = {66, 57};
//+
Line(47) = {65, 56};
//+
Line(48) = {64, 55};
//+
Line(49) = {63, 54};
//+
Line(50) = {62, 53};
//+
Line(51) = {61, 52};
//+
Line(52) = {69, 68};
//+
Line(53) = {68, 67};
//+
Line(54) = {67, 66};
//+
Line(55) = {66, 65};
//+
Line(56) = {65, 64};
//+
Line(57) = {64, 63};
//+
Line(58) = {63, 62};
//+
Line(59) = {62, 61};
//+
Line(60) = {60, 59};
//+
Line(61) = {59, 58};
//+
Line(62) = {58, 57};
//+
Line(63) = {57, 56};
//+
Line(64) = {56, 55};
//+
Line(65) = {55, 54};
//+
Line(66) = {54, 53};
//+
Line(67) = {53, 52};

    // (C) Define periodic boundary condition
    // Periodic Curve {dest} = {src}
        // (1) Connect right boundary {z=Log(R)} to disc boundary.
Periodic Curve {52} = {17};
Periodic Curve {53} = {18};
Periodic Curve {54} = {19};
Periodic Curve {55} = {20};
Periodic Curve {56} = {21};
Periodic Curve {57} = {22};
Periodic Curve {58} = {23};
Periodic Curve {59} = {24};
        // (2) Connect left boundary {z=Log(R_TR)} to {z=Log(R)}.
Periodic Curve {60} = {52};
Periodic Curve {61} = {53};
Periodic Curve {62} = {54};
Periodic Curve {63} = {55};
Periodic Curve {64} = {56};
Periodic Curve {65} = {57};
Periodic Curve {66} = {58};
Periodic Curve {67} = {59};
    // ----

Coherence; // Remove duplicate points

    // -- Define 'Plane Surface' for regions that
    // will be transfinite

Curve Loop(1) = {9, -40, -5, 38};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, -39, -1, 37};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {12, -38, -8, 42};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {4, 37, -8, -41};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {36, 41, -34, -22};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {34, 42, -32, -23};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {29, -19, -27, -39};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {25, 18, -27, 40};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {43, 60, -44, -52};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {44, 61, -45, -53};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {45, 62, -46, -54};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {46, 63, -47, -55};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {47, 64, -48, -56};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {48, 65, -49, -57};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {49, 66, -50, -58};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {50, 67, -51, -59};
//+
Plane Surface(16) = {16};




    // -- Define all 'Transfinite Curve'
        // Elliptical sign-changing interface

        //+
Transfinite Curve {1, 5, 9} = Ceil(Nint/4) Using Progression GeomProgint;
//+
Transfinite Curve {29, 27, 25} = Ceil(Nint/4) Using Progression GeomProgint;
//+
Transfinite Curve {32, 34, 36} = Ceil(Nint/4) Using Progression GeomProgint;
//+
Transfinite Curve {12, 8, 4} = Ceil(Nint/4) Using Progression GeomProgint;
Transfinite Curve {18, 19, 40, 39, 38, 37, 41, 42, 22, 23} = Nmu Using Progression 1;
        // Euler region - horizontal
Transfinite Curve {43, 44, 45, 46, 47, 48, 49, 50, 51} = Nz Using Progression 1;
        // Euler region - vertical
        // {z=Log(R)} vertical is Transfinite. Note that the
        // the number of nodes state here
        // will be overriden by the periodicity condition
Transfinite Curve {52, 53, 54, 55, 56, 57, 58, 59} = 5 Using Progression 1;




    // -- Define all 'Transfinite Surface'

    //+
Transfinite Surface {8} = {36, 37, 74, 77};
//+
Transfinite Surface {7} = {37, 38, 71, 74};
//+
Transfinite Surface {1} = {77, 74, 73, 76};
//+
Transfinite Surface {2} = {74, 71, 70, 73};
//+
Transfinite Surface {5} = {40, 41, 75, 72};
//+
Transfinite Surface {6} = {41, 42, 78, 75};
//+
Transfinite Surface {3} = {75, 78, 76, 73};
//+
Transfinite Surface {4} = {72, 75, 73, 70};
//+
Transfinite Surface {9} = {60, 59, 68, 69};
//+
Transfinite Surface {10} = {59, 58, 67, 68};
//+
Transfinite Surface {11} = {58, 57, 66, 67};
//+
Transfinite Surface {12} = {57, 56, 65, 66};
//+
Transfinite Surface {13} = {56, 55, 64, 65};
//+
Transfinite Surface {14} = {55, 54, 63, 64};
//+
Transfinite Surface {15} = {54, 53, 62, 63};
//+
Transfinite Surface {16} = {53, 52, 61, 62};

    // -- Define 'Plane Surface' for regions
    // that will be unstructured

//+
Curve Loop(17) = {14, 15, 16, 13};
//+
Curve Loop(18) = {25, -17, -24, 32, 12, 9};
//+
Plane Surface(17) = {17, 18};
//+
Curve Loop(19) = {1, 29, 20, 21, 36, 4};
//+
Plane Surface(18) = {19};

    // -- Define Physical Entitites
Physical Surface("omega-m", 68) = {18, 7, 2, 4, 5};
Physical Surface("omega-d", 69) = {1, 8, 6, 3, 17};
Physical Surface("corner-omega-d-pml", 70) = {9, 10, 16, 15};
Physical Surface("corner-omega-m-pml", 71) = {11, 12, 13, 14};
Physical Curve("gamma-d", 72) = {14, 15, 13,16};
Physical Curve("corner-bnd-bot", 73) = {51};
Physical Curve("corner-bnd-right", 74) = {59, 58, 57, 56, 55, 54, 53, 52};

    // -- Mesh parameters
        // Characteristic length based
        // on (approximate) length of transfinite region
length_transfinite = 2*Pi*Sqrt(a_m*b_m);    // (lower bound)
length_transfinite = length_transfinite;
lc = length_transfinite/(Nint);
Printf("Computed characteristic length: %g",lc);
Mesh.CharacteristicLengthMin = lc*CharLengthMin_adim;
Mesh.CharacteristicLengthMax = lc*CharLengthMax_adim;

    // -- Meshing
//Mesh.FlexibleTransfinite = 1;
//Mesh.Algorithm = 6;
//Mesh.RecombinationAlgorithm = 3;
//Mesh.ElementOrder = 2;
If (GenerateQuadMesh>0)
    Include 'Macro-Generate-Quad.geo';
Else
    Mesh 2;
EndIf
