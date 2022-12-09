/*********************************************************************
 * Ellipse perturbed by one corner on the major axis.
 *
 *  - This gmsh file has been partially assembled using the GUI.
 *  - The Euler region is on the disk of radius R_PML, whith R_PML < R.
 *  - Unstructured mesh, structured at the sign-changing interface
 *  - Corner at (x_c,0), symmetric w.r.t. x-axis
 *  - PML region discretized in Euler coordinates (z,theta)
 *  - Recombination algorithm can be used to produce quad mesh
 */

 /* Possible improvements:
 *   - Union of surfaces around sign-changing interface, to avoid the
 *   need for Nint_corner and only rely on Nint.
 *   - Implemenent mesh coarsening using background fields.
 */

// BooleanFragments yields a different line numbering in gmsh 4.9.4 and
// newer, which breaks the geometry.
Warning("This .geo file fails with gmsh 4.9.4 and newer");

SetFactory("OpenCASCADE");
    // ---- Input parameters

    // Ellipse gamma-m semi-axes
If (!Exists(a_m)) a_m=2.5; EndIf
If (!Exists(b_m)) b_m=1; EndIf
    // Ellipse gamma-d semi-axes
If (!Exists(a_d)) a_d=1.3*a_m; EndIf
If (!Exists(b_d)) b_d=Sqrt(a_d^2-Abs(a_m^2-b_m^2)); EndIf
    // Left corner
        // angle in (0,Pi) (rad)
If (!Exists(phi1)) phi1=Pi/2; EndIf
        // R_PML/R radius of PML region
If (!Exists(R1_PML)) R1_PML = 0.7; EndIf
        // R_TR/R radius of truncated region
If (!Exists(R1_TR)) R1_TR = 0.05; EndIf
    // Corner region offset
If (!Exists(x_ofst)) x_ofst = -a_d; EndIf
If (!Exists(y_ofst)) y_ofst = -b_d/2; EndIf
    // Mesh parameters (Structured interface)
        // outer ellipse
If (!Exists(a_o)) a_o=(1-0.1)*a_m+0.1*a_d; EndIf
If (!Exists(b_o)) b_o=Sqrt(Abs(a_o^2-Abs(a_m^2-b_m^2))); EndIf
        // inner ellipse
If (!Exists(a_i)) a_i=(1-2e-2)*a_m; EndIf
If (!Exists(b_i)) b_i = Sqrt(Abs(a_i^2-Abs(a_m^2-b_m^2))); EndIf
        // Sign-changing interface
        // No. of nodes on outside Euler region
If (!Exists(Nint)) Nint=40+1; EndIf
        // No. of nodes in corner disk but outside Euler region
If (!Exists(Nint_corner)) Nint_corner=2; EndIf
        // Geometric progression of element size along interface
If (!Exists(GeomProgint)) GeomProgint=0.5; EndIf
        // No. of nodes in normal direction
If (!Exists(Nmu)) Nmu = 1; EndIf
        // No. of element in z direction
If (!Exists(Nz)) Nz = 5; EndIf
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
p=newl; Ellipse(p) = {0, 0, 0, a_i, b_i, 0, 2*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a_m, b_m, 0, 2*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a_o, b_o, 0, 2*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a_d, b_d, 0, 2*Pi};

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
p = newp; Point(p) = {x_c-R1_PML*R,y_c,0}; // point associated with theta=pi
p = newl; Circle(p) = {x_c, y_c, 0, R, 0, 2*Pi};
p = newl; Circle(p) = {x_c, y_c, 0, R1_PML*R, 0, 2*Pi};
p = newp; Point(p) = {-a_d,0,0};
If (Abs(x_c-R1_PML*R)>a_d)
    Error("Corner circle collides with Gamma-D: reduce R1_PML");
EndIf
    // Junction points with gamma-i, gamma-m, gamma-o
BooleanFragments{ Curve{3}; Curve{2}; Curve{1}; Delete; }{ Curve{5}; Delete; }
Delete {
  Curve{8}; Curve{11}; Curve{14};  Curve{19}; Curve{22}; Curve{16}; Point{19};
}
Line(22) = {5, 11};
Line(23) = {5, 14};
Line(24) = {5, 17};
Line(25) = {5, 18};
Line(26) = {5, 15};
Line(27) = {5, 12};
BooleanFragments{ Curve{22}; Curve{23}; Curve{24}; Curve{25}; Curve{26}; Curve{27}; Delete; }{ Curve{6}; Delete; }
Delete {
  Curve{37};
}
Circle(41) = {20, 19, 6};
Circle(42) = {6, 19, 25};
Delete {
  Curve{22}; Curve{24}; Curve{26}; Curve{32}; Curve{30}; Curve{28};
}
    // Horizontal line {y=0}
Line(43) = {9, 6};
Line(44) = {26, 16};
Line(45) = {16, 13};
Line(46) = {13, 10};
Line(47) = {10, 4};

    // Compute angles of each junction point
xj[] = Point{20}; // gamma-o
phi_o_1 = 2*Atan((xj[1]-y_c)/(xj[0]-x_c));
Printf("[Left corner] Angle outer: %g deg",phi_o_1*(180/Pi));
xj[] = Point{21}; // gamma-m
phi_m_1 = 2*Atan((xj[1]-y_c)/(xj[0]-x_c));
Printf("[Left corner] Angle m: %g deg",phi_m_1*(180/Pi));
xj[] = Point{22}; // gamma-i
phi_i_1 = 2*Atan((xj[1]-y_c)/(xj[0]-x_c));
Printf("[Left corner] Angle inner: %g deg",phi_i_1*(180/Pi));

    // ---- Left corner geometry in Euler coordinates
    // (A) Define points
    // Abscissa
x_input[] = {Log(R*R1_TR), Log(R*R1_PML)};
    // Angle
y_input[] = {-Pi, -phi_o_1/2, -phi1/2,-phi_i_1/2,0,phi_i_1/2,phi1/2,phi_o_1/2,Pi};

For i In {0:(#x_input[]-1)}
    x_input[i] = x_input[i] + x_ofst;
EndFor

For i In {0:(#y_input[]-1)}
    y_input[i] = y_input[i] + y_ofst;
EndFor

For i In {0:(#x_input[]-1)}
    For j In {0:(#y_input[]-1)}
        p = newp; Point(p) = {x_input[i],y_input[j],0};
    EndFor
EndFor

    // (B) Define lines

//+
Line(48) = {36, 37};
//+
Line(49) = {37, 38};
//+
Line(50) = {38, 39};
//+
Line(51) = {39, 40};
//+
Line(52) = {40, 41};
//+
Line(53) = {41, 42};
//+
Line(54) = {42, 43};
//+
Line(55) = {43, 44};
//+
Line(56) = {44, 35};
//+
Line(57) = {43, 34};
//+
Line(58) = {42, 33};
//+
Line(59) = {41, 32};
//+
Line(60) = {40, 31};
//+
Line(61) = {39, 30};
//+
Line(62) = {38, 29};
//+
Line(63) = {37, 28};
//+
Line(64) = {36, 27};
//+
Line(65) = {27, 28};
//+
Line(66) = {28, 29};
//+
Line(67) = {29, 30};
//+
Line(68) = {30, 31};
//+
Line(69) = {31, 32};
//+
Line(70) = {32, 33};
//+
Line(71) = {33, 34};
//+
Line(72) = {34, 35};

    // (C) Define periodic boundary condition
    // Periodic Curve {dest} = {src}
        // (1) Connect right boundary {z=Log(R)} to disc boundary.
Periodic Curve {48} = {42};
Periodic Curve {49} = {38};
Periodic Curve {50} = {39};
Periodic Curve {51} = {40};
Periodic Curve {52} = {34};
Periodic Curve {53} = {35};
Periodic Curve {54} = {36};
Periodic Curve {55} = {41};
        // (2) Connect left boundary {z=Log(R_TR)} to {z=Log(R)}.
Periodic Curve {72} = {55};
Periodic Curve {71} = {54};
Periodic Curve {70} = {53};
Periodic Curve {69} = {52};
Periodic Curve {68} = {51};
Periodic Curve {67} = {50};
Periodic Curve {66} = {49};
Periodic Curve {65} = {48};

    // ---

    // -- Define 'Plane Surface' for regions that
    // will be transfinite
//+
Curve Loop(1) = {36, 23, -18, -25};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {18, -7, -46, 10};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {38, 31, -20, -33};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {12, 46, -9, 20};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {39, 29, -21, -31};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {21, 15, 45, -12};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {45, 10, -17, -13};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {56, -72, -57, 55};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {54, 57, -71, -58};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {53, 58, -70, -59};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {52, 59, -69, -60};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {51, 60, -68, -61};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {50, 61, -67, -62};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {49, 62, -66, -63};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {48, 63, -65, -64};
//+
Plane Surface(15) = {15};

Curve Loop(16) = {25, -17, -27, 35};
//+
Plane Surface(16) = {16};

        // -- Define all 'Transfinite Curve'
        // Elliptical sign-changing interface
Transfinite Curve {7, 10, 13,15, 12, 9} = Nint/2 Using Bump GeomProgint;
Transfinite Curve {23, 25, 27, 33, 31, 29} = Nint_corner Using Progression 1;
Transfinite Curve {36, 35, 18, 17, 39, 38, 20, 21, 45, 46} = Nmu Using Progression 1;
        // Euler region - horizontal
Transfinite Curve {56, 57, 58, 59, 60, 61, 62, 63, 64} = Nz Using Progression 1;
        // Euler region - vertical
        // {z=Log(R)} vertical is Transfinite. Note that the
        // the number of nodes state here
        // will be overriden by the eriodicity condition
Transfinite Curve {55, 54, 53, 52, 51, 50, 49, 48, 72, 71, 70, 69, 68, 67, 66, 65} = 5 Using Progression 1;

        // -- Define all 'Transfinite Surface'
//+
Transfinite Surface {1} = {20, 21, 14, 11};
//+
Transfinite Surface {16} = {21, 22, 17, 14};
//+
Transfinite Surface {7} = {14, 17, 16, 13};
//+
Transfinite Surface {2} = {11, 14, 13, 10};
//+
Transfinite Surface {5} = {23, 24, 15, 18};
//+
Transfinite Surface {6} = {18, 15, 13, 16};
//+
Transfinite Surface {3} = {24, 25, 12, 15};
//+
Transfinite Surface {4} = {15, 12, 10, 13};
//+
Curve Loop(17) = {4};
//+
//+
Transfinite Surface {8} = {34, 43, 44, 35};
//+
Transfinite Surface {9} = {33, 42, 43, 34};
//+
Transfinite Surface {10} = {32, 41, 42, 33};
//+
Transfinite Surface {11} = {31, 40, 41, 32};
//+
Transfinite Surface {12} = {30, 39, 40, 31};
//+
Transfinite Surface {13} = {29, 38, 39, 30};
//+
Transfinite Surface {14} = {28, 37, 38, 29};
//+
Transfinite Surface {15} = {27, 36, 37, 28};

        // -- Define 'Plane Surface' for regions
        // that will be unstructured

// Ellipse Gamma_D should be two ellipse arcs
BooleanFragments{ Curve{4}; Delete; }{ Curve{43}; Curve{47}; }

Curve Loop(18) = {73, 43, -41, 23, -7, 47};
Plane Surface(17) = {18};
Curve Loop(19) = {43, 42, 33, 9, 47, -74};
Curve Loop(20) = {42, 33, 9, 47, -74, 43};
Plane Surface(18) = {20};
Curve Loop(21) = {13, -27, -34, 44};
Plane Surface(19) = {21};
Curve Loop(22) = {44, -15, -29, 40};
Plane Surface(20) = {22};

        // -- Define Physical Entitites
//+
Physical Curve("gamma-d", 75) = {73, 74};
//+
Physical Surface("omega-d", 76) = {17, 2, 1, 18, 4, 3};
//+
Physical Surface("omega-m", 77) = {7, 19, 20, 6, 5, 16};
//+
Physical Curve("corner-bnd-bot", 78) = {64};
//+
Physical Curve("corner-bnd-right", 79) = {48, 49, 50, 51, 52, 53, 54, 55};
//+
Physical Surface("corner-omega-d-pml", 80) = {15, 14, 9, 8};
//+
Physical Surface("corner-omega-m-pml", 81) = {10, 11, 12, 13};

        // -- Mesh parameters
        // Characteristic length based
        // on (approximate) length of transfinite region
length_transfinite = 2*Pi*Sqrt(a_m*b_m);    // (lower bound)
length_transfinite = length_transfinite - R*R1_PML;
lc = length_transfinite/(Nint);
Printf("Computed characteristic length: %g",lc);
Mesh.CharacteristicLengthMin = lc*CharLengthMin_adim;
Mesh.CharacteristicLengthMax = lc*CharLengthMax_adim;

        // -- Meshing

//Mesh.Algorithm = 6;
//Mesh.RecombinationAlgorithm = 2;
//Mesh.ElementOrder = 2;
If (GenerateQuadMesh>0)
    Include 'Macro-Generate-Quad.geo';
Else
    Mesh 2;
EndIf

/*
// Original meshing using refinement.
// Issue: this does not preserve Nmu in boundary layer region.
If (GenerateQuadMesh)
    // This methodology can leave triangle around the corner.
        // Mesh with Frontal-Delaunay (Default)
    Mesh 2;
        // Recombine 2D with Blossom (Default)
    RecombineMesh;
        // Refine by splitting
    RefineMesh;
        // Recombine 2D with Blossom (Default)
    RecombineMesh;
    Printf("[Statistics] Number of remaining triangles: %g",Mesh.NbTriangles);
    If (Mesh.NbTriangles>0)
        Error("There are %g remaining triangles in the mesh. Adjust the distribution of nodes along transfinite sign-changing interface (typically Nint_corner).",Mesh.NbTriangles);
    EndIf
EndIf

*/
