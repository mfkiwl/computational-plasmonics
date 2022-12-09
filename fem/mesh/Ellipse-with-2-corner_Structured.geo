/*********************************************************************
 * Ellipse perturbed by two corners, one on the major axis and one on
 * the minor axis. For both corners, the Euler region is on the disk of
 * radius R.
 * This gmsh file has been partially assembled using the GUI.
 *  - Unstructured mesh, structured at the sign-changing interface.
 *  - Corner at (x_c,0), symmetric w.r.t. x-axis
 *  - Corner at (0,y_c), symmetric w.r.t. y-axis
 *  - Corner disks discretized in Euler coordinates (z,theta)
 *  - Recombination algorithm can be used to produce quad mesh.
 */

 // TODO: Adjust element distribution on sign-changing interface (too dense on
 // top left boundary.)

// BooleanFragments yields a different line numbering in gmsh 4.9.4 and
// newer, which breaks the geometry.
Warning("This .geo file fails with gmsh 4.9.4 and newer");

SetFactory("OpenCASCADE");
    // ---- Input parameters

    // Ellipse gamma-m semi-axes
If (!Exists(a_m)) a_m=2.5; EndIf
If (!Exists(b_m)) b_m=2.45; EndIf
c = Abs(a_m^2-b_m^2);
    // Ellipse gamma-d semi-axes
If (!Exists(a_d)) a_d=1.7*a_m; EndIf
If (!Exists(b_d)) b_d=Sqrt(a_d^2-c); EndIf
    // Left corner
If (!Exists(phi1)) phi1=0.75*Pi; EndIf // angle in (0,Pi) (rad)
If (!Exists(R1_TR)) R1_TR = 0.05; EndIf // R_TR/R radius of truncated region
    // Top corner
If (!Exists(phi2)) phi2=0.75*Pi; EndIf // angle in (0,Pi) (rad)
If (!Exists(R2_TR)) R2_TR = 0.05; EndIf // R_TR/R radius of truncated region
    // Corner region offset
If (!Exists(x_ofst_1)) x_ofst_1 = -a_d-0.5; EndIf
If (!Exists(y_ofst_1)) y_ofst_1 = -b_d; EndIf
If (!Exists(x_ofst_2)) x_ofst_2 = -a_d-0.5; EndIf
If (!Exists(y_ofst_2)) y_ofst_2 = y_ofst_1+2*Pi+0.1; EndIf

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

    // ---- Points associated with corner on minor axis (C^1 junction)
    // Junction point on gamma-m
    // (slope tan(psi2/2) and yj>0)

psi2 = Pi - phi2;
x_j = Sqrt(Tan(psi2/2)^2+(b_m/a_m)^2);
x_j = -a_m * Tan(psi2/2) / x_j;
y_j = b_m*Sqrt(1-(x_j/a_m)^2);
    // corner and radius
x_c = 0;
y_c = y_j - x_j*Tan(psi2/2);
p = newp; Point(p) = {x_c,y_c,0};
R = Sqrt((x_c-x_j)^2+(y_c-y_j)^2);
p = newl; Circle(p) = {x_c, y_c, 0, R, -Pi/2, Pi/2};
p = newl; Circle(p) = {x_c, y_c, 0, R, Pi/2, 3*Pi/2};
If (Abs(y_c+R)>b_d)
    Error("Corner circle collides with Gamma-D: increase phi2");
EndIf
R2 = R;
    // Junction points
BooleanFragments{ Curve{38}; Curve{37}; Delete; }{ Curve{29}; Curve{1}; Curve{5}; Curve{9}; Curve{27}; Curve{25}; Delete; }
//+
Recursive Delete {
  Curve{55}; Curve{53}; Curve{45}; Curve{48}; Curve{50}; Curve{52};
}

    // Angle of each junction point
xj[] = Point{58}; // gamma-o
phi_o = 2*Atan((xj[0]-x_c)/(y_c-xj[1]));
Printf("[Top corner] Angle outer: %g deg",phi_o*(180/Pi));
xj[] = Point{59}; // gamma-m
phi_m = 2*Atan((xj[0]-x_c)/(y_c-xj[1]));
Printf("[Top corner] Angle m: %g deg",phi_m*(180/Pi));
xj[] = Point{60}; // gamma-i
phi_i = 2*Atan((xj[0]-x_c)/(y_c-xj[1]));
Printf("[Top corner] Angle inner: %g deg",phi_i*(180/Pi));
phi2_o = phi_o; phi2_m = phi_m; phi2_i = phi_i;

    // ---- Lines to subdivise interface boundary layer
Line(57) = {8, 16};
//+
Line(58) = {16, 24};
//+
Line(59) = {9, 17};
//+
Line(60) = {17, 25};

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
Line(61) = {83, 74};
//+
Line(62) = {82, 73};
//+
Line(63) = {81, 72};
//+
Line(64) = {80, 71};
//+
Line(65) = {79, 70};
//+
Line(66) = {78, 69};
//+
Line(67) = {77, 68};
//+
Line(68) = {76, 67};
//+
Line(69) = {75, 66};
//+
Line(70) = {83, 82};
//+
Line(71) = {82, 81};
//+
Line(72) = {81, 80};
//+
Line(73) = {80, 79};
//+
Line(74) = {79, 78};
//+
Line(75) = {78, 77};
//+
Line(76) = {77, 76};
//+
Line(77) = {76, 75};
//+
Line(78) = {74, 73};
//+
Line(79) = {73, 72};
//+
Line(80) = {72, 71};
//+
Line(81) = {71, 70};
//+
Line(82) = {70, 69};
//+
Line(83) = {69, 68};
//+
Line(84) = {68, 67};
//+
Line(85) = {67, 66};

    // (C) Define periodic boundary condition
    // Periodic Curve {dest} = {src}
        // (1) Connect right boundary {z=Log(R)} to disc boundary.
Periodic Curve {70} = {24};
Periodic Curve {71} = {23};
Periodic Curve {72} = {22};
Periodic Curve {73} = {21};
Periodic Curve {74} = {20};
Periodic Curve {75} = {19};
Periodic Curve {76} = {18};
Periodic Curve {77} = {17};
        // (2) Connect left boundary {z=Log(R_TR)} to {z=Log(R)}.
Periodic Curve {78} = {70};
Periodic Curve {79} = {71};
Periodic Curve {80} = {72};
Periodic Curve {81} = {73};
Periodic Curve {82} = {74};
Periodic Curve {83} = {75};
Periodic Curve {84} = {76};
Periodic Curve {85} = {77};
    // ----

    // ---- Top corner geometry in Euler coordinates
    // (A) Define points
    // Abscissa
x_input[] = {Log(R2*R2_TR), Log(R2)};
    // Angle
y_input[] = {-Pi, -phi2_o/2, -phi2_m/2,-phi2_i/2,0,phi2_i/2,phi2_m/2,phi2_o/2,Pi};

For i In {0:(#x_input[]-1)}
    x_input[i] = x_input[i] + x_ofst_2;
EndFor

For i In {0:(#y_input[]-1)}
    y_input[i] = y_input[i] + y_ofst_2;
EndFor

For i In {0:(#x_input[]-1)}
    For j In {0:(#y_input[]-1)}
        p = newp; Point(p) = {x_input[i],y_input[j],0};
    EndFor
EndFor

    // (B) Define lines
//+
Line(86) = {100, 91};
//+
Line(87) = {99, 90};
//+
Line(88) = {98, 89};
//+
Line(89) = {97, 88};
//+
Line(90) = {96, 87};
//+
Line(91) = {95, 86};
//+
Line(92) = {94, 85};
//+
Line(93) = {93, 84};
//+
Line(94) = {101, 92};
//+
Line(95) = {101, 100};
//+
Line(96) = {100, 99};
//+
Line(97) = {99, 98};
//+
Line(98) = {98, 97};
//+
Line(99) = {97, 96};
//+
Line(100) = {96, 95};
//+
Line(101) = {95, 94};
//+
Line(102) = {94, 93};
//+
Line(103) = {92, 91};
//+
Line(104) = {91, 90};
//+
Line(105) = {90, 89};
//+
Line(106) = {89, 88};
//+
Line(107) = {88, 87};
//+
Line(108) = {87, 86};
//+
Line(109) = {86, 85};
//+
Line(110) = {85, 84};

    // (C) Define periodic boundary condition
    // Periodic Curve {dest} = {src}
        // (1) Connect right boundary {z=Log(R)} to disc boundary.
Periodic Curve {95} = {41};
Periodic Curve {96} = {42};
Periodic Curve {97} = {43};
Periodic Curve {98} = {44};
Periodic Curve {99} = {37};
Periodic Curve {100} = {38};
Periodic Curve {101} = {39};
Periodic Curve {102} = {40};
        // (2) Connect left boundary {z=Log(R_TR)} to {z=Log(R)}.
Periodic Curve {103} = {95};
Periodic Curve {104} = {96};
Periodic Curve {105} = {97};
Periodic Curve {106} = {98};
Periodic Curve {107} = {99};
Periodic Curve {108} = {100};
Periodic Curve {109} = {101};
Periodic Curve {110} = {102};

    // ----

Coherence; // Remove duplicate points

    // -- Define 'Plane Surface' for regions that
    // will be transfinite

//+
Curve Loop(1) = {49, -42, -51, -60};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {43, -47, 59, 49};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {39, 56, 18, -54};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {38, 54, 19, -46};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {22, 34, -57, -36};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {23, 32, -58, -34};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {8, 60, -12, -58};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {4, 59, -8, -57};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {94, 103, -86, -95};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {86, 104, -87, -96};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {87, 105, -88, -97};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {88, 106, -89, -98};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {89, 107, -90, -99};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {90, 108, -91, -100};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {91, 109, -92, -101};
//+
Plane Surface(15) = {15};
//+
Curve Loop(16) = {92, 110, -93, -102};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {61, 78, -62, -70};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {62, 79, -63, -71};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {63, 80, -64, -72};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {64, 81, -65, -73};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {65, 82, -66, -74};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {66, 83, -67, -75};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {67, 84, -68, -76};
//+
Plane Surface(23) = {23};
//+
Curve Loop(24) = {68, 85, -69, -77};
//+
Plane Surface(24) = {24};

    // -- Define all 'Transfinite Curve'
        // Elliptical sign-changing interface
//+
Transfinite Curve {47, 49, 51} = Ceil(Nint/4) Using Progression GeomProgint;
//+
Transfinite Curve {46, 54, 56} = Ceil(Nint/4) Using Progression 1;
//+
Transfinite Curve {4, 8, 12} = Ceil(Nint/4) Using Progression 1/GeomProgint;
//+
Transfinite Curve {32, 34, 36} = Ceil(Nint/4) Using Progression 1;
Transfinite Curve {39, 18, 19, 38, 42, 43, 60, 59, 57, 58, 23, 22} = Nmu Using Progression 1;
        // Euler region - horizontal
Transfinite Curve {94, 86, 87, 88, 89, 90, 91, 92, 93, 61, 62, 63, 64, 65, 66, 67, 68, 69} = Nz Using Progression 1;

        // Euler region - vertical
        // {z=Log(R)} vertical is Transfinite. Note that the
        // the number of nodes state here
        // will be overriden by the periodicity condition
Transfinite Curve {95, 96, 97, 98, 99, 100, 101, 102, 70, 71, 72, 73, 74, 75, 76, 77} = 5 Using Progression 1;

    // -- Define all 'Transfinite Surface'
//+
Transfinite Surface {3} = {36, 37, 55, 56};
//+
Transfinite Surface {4} = {37, 38, 54, 55};
//+
Transfinite Surface {1} = {59, 105, 107, 58};
//+
Transfinite Surface {2} = {60, 103, 105, 59};
//+
Transfinite Surface {8} = {103, 102, 104, 105};
//+
Transfinite Surface {7} = {104, 106, 107, 105};
//+
Transfinite Surface {5} = {40, 41, 104, 102};
//+
Transfinite Surface {6} = {41, 42, 106, 104};
//+
Transfinite Surface {9} = {92, 91, 100, 101};
//+
Transfinite Surface {10} = {91, 90, 99, 100};
//+
Transfinite Surface {11} = {90, 89, 98, 99};
//+
Transfinite Surface {12} = {89, 88, 97, 98};
//+
Transfinite Surface {13} = {88, 87, 96, 97};
//+
Transfinite Surface {14} = {87, 86, 95, 96};
//+
Transfinite Surface {15} = {86, 85, 94, 95};
//+
Transfinite Surface {16} = {85, 84, 93, 94};
//+
Transfinite Surface {17} = {74, 73, 82, 83};
//+
Transfinite Surface {18} = {73, 72, 81, 82};
//+
Transfinite Surface {19} = {72, 71, 80, 81};
//+
Transfinite Surface {20} = {71, 70, 79, 80};
//+
Transfinite Surface {21} = {70, 69, 78, 79};
//+
Transfinite Surface {22} = {69, 68, 77, 78};
//+
Transfinite Surface {23} = {68, 67, 76, 77};
//+
Transfinite Surface {24} = {67, 66, 75, 76};


    // -- Define 'Plane Surface' for regions
    // that will be unstructured
//+
Curve Loop(25) = {46, 20, 21, 36, 4, 47, 44, 37};
//+
Plane Surface(25) = {25};
//+
Curve Loop(26) = {41, -51, -12, -32, 24, 17, -56, 40};
//+
Curve Loop(27) = {15, 16, 13, 14};
//+
Curve Loop(28) = {40, 41, -51, -12, -32, 24, 17, -56};
//+
Plane Surface(26) = {27, 28};


    // -- Define Physical Entitites
Physical Curve("gamma-d", 111) = {14, 13, 16, 15};
//+
Physical Surface("omega-d", 112) = {26, 3, 1, 7, 6};
//+
Physical Surface("omega-m", 113) = {4, 2, 8, 5, 25};
//+
Physical Surface("corner-1-omega-d-pml", 114) = {17, 18, 24, 23};
//+
Physical Surface("corner-2-omega-d-pml", 115) = {16, 15, 9, 10};
//+
Physical Surface("corner-2-omega-m-pml", 116) = {11, 12, 13, 14};
//+
Physical Surface("corner-1-omega-m-pml", 117) = {19, 20, 21, 22};
//+
Physical Curve("corner-2-bnd-bot", 118) = {93};
//+
Physical Curve("corner-1-bnd-bot", 119) = {69};
//+
Physical Curve("corner-1-bnd-right", 120) = {77, 76, 75, 74, 73, 72, 71, 70};
//+
Physical Curve("corner-2-bnd-right", 121) = {102, 101, 100, 99, 98, 97, 96, 95};

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
//Mesh.Algorithm = 8;
//Mesh.RecombinationAlgorithm = 3;
Mesh.ElementOrder = 2;
If (GenerateQuadMesh>0)
    Include 'Macro-Generate-Quad.geo';
Else
    Mesh 2;
EndIf
