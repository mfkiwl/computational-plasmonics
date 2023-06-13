/*********************************************************************
 * Unstructured mesh of an ellipse.
 *  - No symmetry is enforced.
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
If (!Exists(Nint)) Nint=20; EndIf
        // Geometric progression of element size along interface
If (!Exists(GeomProgint)) GeomProgint=0.4; EndIf
    // Adimensional characteristic length
If (!Exists(CharLengthMin_adim)) CharLengthMin_adim=1; EndIf
If (!Exists(CharLengthMax_adim)) CharLengthMax_adim=4; EndIf

    // Recombination algorithm are very sensitive to mesh size.
    // If recombination fails or is only partial (i.e. there are
    // triangles remaining) then tweak mesh size at interface and/or
    // characteristic mesh size.
    // See 'Unstructured quadrangular meshes' gmsh tutorial.
If (!Exists(GenerateQuadMesh)) GenerateQuadMesh=1; EndIf

    // ---- Ellipses
p=newl; Ellipse(p) = {0, 0, 0, a_d, b_d, 0, 2*Pi};
p=newl; Ellipse(p) = {0, 0, 0, a_m, b_m, 0, 2*Pi};
p=newp; Point(p) = {a_m,0,0,0};
p=newp; Point(p) = {-a_m,0,0,0};
    // cut ellipse omega-m in two (for orthoradial refinment)
//+
Line(3) = {4, 2};
//+
BooleanFragments{ Curve{2}; Delete; }{ Curve{3}; Delete; }
//+
Recursive Delete {
  Curve{3};
}

    // -- Define all 'Transfinite Curve'
    // Elliptical sign-changing interface
Transfinite Curve {4,5} = Floor(Nint/2) Using Bump GeomProgint;

    // -- Define 'Plane Surface' for regions
    // that will be unstructured
//+
Curve Loop(1) = {4, 5};
//+
Curve Loop(2) = {4, 5};
//+
Plane Surface(1) = {2};
//+
Curve Loop(3) = {1};
//+
Curve Loop(4) = {4, 5};
//+
Plane Surface(2) = {3, 4};
    // -- Define Physical Entitites
//+
Physical Curve("gamma-d", 6) = {1};
//+
Physical Surface("omega-d", 7) = {2};
//+
Physical Surface("omega-m", 8) = {1};
//+
Physical Curve("gamma-m", 7) = {4,5};

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
