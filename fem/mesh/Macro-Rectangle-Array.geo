/* Array of transfinite rectangles.
 * 
 * This geometry can be imported in a
 * parent geometry using the 'Import' statement.
 * The input parameters are listed below along their
 * default values.
 * 
 * To debug, display the mesh in gmsh with the labels.

Positive orientation:
---------------------
y
^
|
|
+---> x

Numbering convention:
---------------------
                  .
                  .
                  .
      p(k+N_x)     p(k+N_x+1)    
      +-----------+-----------+
      |           |           |
      |   s(k)    |  s(k+1)   |
. . . |           ^ L_k+1,    | . . .
      |           |   k+N_x+1 |
      |   L_k,k+1 |           |
      +----->-----+-----------+
      p(k)        p(k+1)      p(k+2)
                  .
                  .
                  .
                  
- Points: linear indexing starting from
    bottom left.
- Edges: L_k,p is the *positively oriented*
    edge from k to p.
- Rectangle: s(k) is the rectangle (2D) whose
    bottom left point is p(k). This numbering
    is non-continuous.
*/

SetFactory("OpenCASCADE");

    // ---- Input parameters

        // Domain name
        // string
If (!Exists(name)) name = "rectangle"; EndIf

        // x coordinates of points
        // must be a strictly increasing sequence
If (!Exists(x_input)) x_input[] = { 0, 1, 2}; EndIf

        // y coordinates of points
        // must be a strictly increasing sequence
If (!Exists(y_input)) y_input[] = { 0, 1, 2, 3}; EndIf

        // Number of nodes in x-direction
        // Length N_x-1
If (!Exists(N_mesh_x))
    For i In {0:(#x_input[]-2)}
        N_mesh_x[i] = 5;
    EndFor
EndIf
        // Number of nodes in y-direction
        // Length N_y-1
If (!Exists(N_mesh_y))
    For i In {0:(#y_input[]-2)}
        N_mesh_y[i] = 5;
    EndFor
EndIf
        // Define physical entities
        // If a parent geometry already defines physical entities in
        // the rectangle, set this setting to '0' to avoid overlapping
        // physical entities. (Overlapping physical entities can cause
        // some FEM packages to fail.)

If (!Exists(Bool_define_line_physical_entities))
    Bool_define_line_physical_entities = 1;
EndIf

If (!Exists(Bool_define_surface_physical_entities))
    Bool_define_surface_physical_entities = 1;
EndIf

    // ---- Validate input arguments
N_x = #x_input[];
N_y = #y_input[];

        // Check sizes
If ( ((#N_mesh_x[])!=(N_x-1)) || ((#N_mesh_y[])!=(N_y-1)))
    Error("Input coordinates and mesh size vectors do not match.");
EndIf
    // Check monotony
For i In {0:(N_x-2)}
    If (x_input[i]>x_input[i+1])
        Error("Input x-coordinates vector is not strictly increasing.");
    EndIf
EndFor

For j In {0:(N_y-2)}
    If (y_input[j]>y_input[j+1])
        Error("Input y-coordinates vector is not strictly increasing.");
    EndIf
EndFor

    // ---- Build geometry

Macro LinearIndex
    // Linear index of point #(i,j)
    k = i+j*N_x;
Return

    // Build all points
N_p = N_x*N_y; // number of points
For j In {0:(N_y-1)}
    For i In {0:(N_x-1)}
        Call LinearIndex;
        x[k] = x_input[i]; y[k] = y_input[j];
        p[k] = newp; Point(p[k]) = {x[k],y[k],0};
        //Printf("Point(%g) = (%g,%g)",k,i,j);
    EndFor
EndFor

    // Build all lines
For j In {0:(N_y-1)}
    For i In {0:(N_x-1)}
        Call LinearIndex;
        If (i<(N_x-1))
            // horizontal positively oriented edge
            l~{k}~{k+1} = newl; Line(l~{k}~{k+1}) = {p[k],p[k+1]};
            Transfinite Curve {l~{k}~{k+1}} = N_mesh_x[i] Using Progression 1;
        EndIf
        If (j<(N_y-1))
            // vertical positively oriented edge
            l~{k}~{k+N_x} = newl; Line(l~{k}~{k+N_x}) = {p[k],p[k+N_x]};
            Transfinite Curve {l~{k}~{k+N_x}} = N_mesh_y[j] Using Progression 1;
        EndIf
    EndFor
EndFor


    // Build all rectangles
        // Loop over bottom-left point of the rectangle
For j In {0:(N_y-2)}
    For i In {0:(N_x-2)}
    Call LinearIndex;
    l_Rect = newl;
    Line Loop(l_Rect) = {l~{k}~{k+1},l~{k+1}~{k+1+N_x},-l~{k+N_x}~{k+N_x+1},-l~{k}~{k+N_x}};
    s~{k} = news;
    Plane Surface(s~{k}) = {l_Rect};
    EndFor
EndFor

    // ---- Physical entities
If (Bool_define_line_physical_entities == 1)
    Physical Line(StrCat(name,"-bnd-right"))={};
    i=N_x-1;
    For j In {0:(N_y-2)}
        Call LinearIndex;
        Physical Line(StrCat(name,"-bnd-right")) += {l~{k}~{k+N_x}};
    EndFor

    Physical Line(StrCat(name,"-bnd-left"))={};
    i=0;
    For j In {0:(N_y-2)}
        Call LinearIndex;
        Physical Line(StrCat(name,"-bnd-left")) += {l~{k}~{k+N_x}};
    EndFor

    Physical Line(StrCat(name,"-bnd-bot"))={};
    j=0;
    For i In {0:(N_x-2)}
        Call LinearIndex;
        Physical Line(StrCat(name,"-bnd-bot")) += {l~{k}~{k+1}};
    EndFor

    Physical Line(StrCat(name,"-bnd-top"))={};
    j=N_y-1;
    For i In {0:(N_x-2)}
        Call LinearIndex;
        Physical Line(StrCat(name,"-bnd-top")) += {l~{k}~{k+1}};
    EndFor
EndIf

    // Surfaces
idx_surface = 1; // 1-based numbering of physical entities
For j In {0:(N_y-2)}
    For i In {0:(N_x-2)}
    Call LinearIndex;
    If (Bool_define_surface_physical_entities == 1)
        Physical Surface(StrCat(name,Sprintf("-surf-%02g",k))) = {s~{k}};
    EndIf
    Transfinite Surface {s~{k}};
    idx_surface += 1;
    EndFor
EndFor


Macro ConnectSrcToRightBnd
    // Connect right boundary to l_src.
    // l_src: array of lines (length N_y-1)
    
        // Loop over all vertical edge
    For j In {0:(N_y-2)}
        For i In {0:(N_x-1)}
            Call LinearIndex;
                // vertical positively oriented edge
                Periodic Curve {l~{k}~{k+N_x}} = {l_src[j]};
        EndFor
    EndFor    
Return
