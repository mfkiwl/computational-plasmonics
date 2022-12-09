/* Recombine triangles into quads.
*
* The recombination strategy is based on the value of 'GenerateQuadMesh':
*   1: Mesh, Recombine (recombination includes a subdivision step).
*   2: Mesh, Subdivide.
* See tutorial 't11.geo' and documentation.
*
* The result of the recombination (in particular the remaining number of
* triangles) can be affected by the following parameters, to be set by the
* calling script:
*   - Mesh.Algorithm (Default: Frontal-Delaunay)
*       Meshing algorithm, affects (1,2).
*       '8: Frontal-Delaunay for Quads' and '9: Packing of Parallelograms' are
*       designed for producing quads.
*   - Mesh.RecombinationAlgorithm (Default: Blossom)
*       Recombination algorithm, affects only (1).
*       '2: simple full-quad' and '3: blossom full-quad' are designed for
*       producing quads.
*   - Mesh.SubdivisionAlgorithm (Default: 0)
*       Algorithm used by subdivision step, affects (1,2).
*   - Mesh.FlexibleTransfinite
*
*  Recombination of triangles into quads is very sensitive to the initial mesh.
*  If some triangles are remaining, tweak initial mesh.
*/

RecombineAll = Mesh.RecombineAll;

Printf("[Macro-Generate-Quad] Mesh.Algorithm = %g",Mesh.Algorithm);
Printf("[Macro-Generate-Quad] Mesh.RecombinationAlgorithm = %g",Mesh.RecombinationAlgorithm);
Printf("[Macro-Generate-Quad] Mesh.SubdivisionAlgorithm = %g",Mesh.SubdivisionAlgorithm);

If (GenerateQuadMesh==1)
    Printf("[Macro-Generate-Quad] (1) Mesh and recombine.");
        // Recombination will be done after *each* surface is meshed.
        // Mesh.RecombineAll=0 leads to a sequential approach where recombination
        // is performed after *all* surfaces have been meshed, which is typically
        // less likely to recombine all triangles.
    Mesh.RecombineAll = 1;
    Mesh 2;
ElseIf (GenerateQuadMesh==2)
    Printf("[Macro-Generate-Quad] (2) Mesh and subdivide.");
        // (a) Mesh
    Mesh.RecombineAll = 0;
    Mesh 2; // Mesh
        // (b) Refine by subdivision (order reduced to 1)
    order = Mesh.ElementOrder;
    SubdivisionAlgorithm = Mesh.SubdivisionAlgorithm;
    Mesh.SubdivisionAlgorithm = 1;
    RefineMesh;
    Mesh.SubdivisionAlgorithm = SubdivisionAlgorithm;
        // (c) Restore original order
    SetOrder order;
EndIf

If (Mesh.NbTriangles>0)
    Warning("[Macro-Generate-Quad] There are %g remaining triangles in the mesh.",Mesh.NbTriangles);
Else
    Printf("[Macro-Generate-Quad] All triangles have been recombined into quads.");
EndIf
    // Restore modified options
Mesh.RecombineAll = RecombineAll;
