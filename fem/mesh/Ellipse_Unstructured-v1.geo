/*********************************************************************
 * Unstructured mesh of an ellipse.
 *  - Major and minor axes added to enforce symmetries.
 *  - Recombination algorithm can be used to produce quad mesh.
 */

SetFactory("OpenCASCADE");
    // ---- Input parameters

    // Ellipse gamma-m semi-axes
If (!Exists(a_m)) a_m=2.5; EndIf
If (!Exists(b_m)) b_m=1; EndIf
    // Ellipse gamma-d semi-axes
If (!Exists(a_d)) a_d=1.3*a_m; EndIf
If (!Exists(b_d)) b_d=Sqrt(a_d^2-Abs(a_m^2-b_m^2)); EndIf
        // Sign-changing interface
        // No. of nodes
If (!Exists(Nint)) Nint=90; EndIf
        // Geometric progression of element size along interface
If (!Exists(GeomProgint)) GeomProgint=1.01; EndIf
    // Adimensional characteristic length
If (!Exists(CharLengthMin_adim)) CharLengthMin_adim=1; EndIf
If (!Exists(CharLengthMax_adim)) CharLengthMax_adim=2; EndIf
If (!Exists(GenerateQuadMesh)) GenerateQuadMesh=1; EndIf

    // ---- Ellipses
p=newl; Ellipse(p) = {0, 0, 0, a_d, b_d, 0, 2*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a_m, b_m, 0, 2*Pi};
    // cut ellipse-m in four
p=newp; Point(p) = {a_d,0,0,0};
p=newp; Point(p) = {-a_d,0,0,0};
p=newp; Point(p) = {0,b_d,0,0};
p=newp; Point(p) = {0,-b_d,0,0};
Line(3) = {4, 1};
Line(4) = {5, 6};
BooleanFragments{ Curve{2}; Delete; }{ Curve{4}; Curve{3}; Delete; }
    // cut ellipse-d in two
BooleanFragments{ Curve{1}; Delete; }{ Curve{6}; Curve{9}; Curve{11}; Curve{10}; Delete; }

    // -- Define all 'Transfinite Curve'
    // Elliptical sign-changing interface
Transfinite Curve {2} = Ceil(Nint/4) Using Progression GeomProgint;
Transfinite Curve {3} = Ceil(Nint/4) Using Progression 1/GeomProgint;
Transfinite Curve {4} = Ceil(Nint/4) Using Progression GeomProgint;
Transfinite Curve {5} = Ceil(Nint/4) Using Progression 1/GeomProgint;

    // -- Define 'Plane Surface' for regions
    // that will be unstructured
//+
Curve Loop(1) = {3, 12, -7};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, 7, 13};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {12, 8, -4};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {13, -5, -8};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {14, 6, -2, 10};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {17, -10, -5, 9};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {11, 4, 9, -16};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {15, 11, -3, -6};
//+
Plane Surface(8) = {8};

    // -- Define Physical Entitites
Physical Curve("gamma-d", 18) = {15, 14, 17, 16};
Physical Surface("omega-d", 19) = {8, 5, 6, 7};
Physical Surface("omega-m", 20) = {1, 2, 4, 3};

    // -- Meshing
        // Characteristic length based
        // on (approximate) length of transfinite region
length_transfinite = 2*Pi*Sqrt(a_m*b_m);    // (lower bound)
lc = length_transfinite/(Nint);
Printf("Computed characteristic length: %g",lc);
Mesh.CharacteristicLengthMin = lc*CharLengthMin_adim;
Mesh.CharacteristicLengthMax = lc*CharLengthMax_adim;

//Mesh.RecombinationAlgorithm = 3;
//Mesh.Algorithm = 9;
//Mesh.ElementOrder = 2;
If (GenerateQuadMesh>0)
    Include 'Macro-Generate-Quad.geo';
Else
    Mesh 2;
EndIf
