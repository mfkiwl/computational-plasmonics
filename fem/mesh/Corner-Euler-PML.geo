/* Corner geometry 
 *      {R_TR<r<R, -pi<theta<pi}
 * in Euler coordinates (z,theta).
 * The computational region is:
 *      {ln(R_TR) < z < ln(R), -pi < theta < pi}
 * The PML region is:
 *      {ln(R_TR) < z < ln(R_PML), -pi < theta < pi}
 * 
 * Physical entities:
 *  omega-m: particle outside PML region
 *  omega-d: dielectric outside PML region
 *  gamma-d: Dirichlet boundary
 *  gamma-top: top boundary theta=+pi
 *  gamma-bot: bottom boundary theta=-pi
 *  omega-m-pml: PML region inside omega-m 
 *  omega-d-pml: PML region inside omega-d 
*/
SetFactory("OpenCASCADE");
    // ---- Input parameters
        // Abscissa (z coordinates)
If (!Exists(R)) R=1; EndIf
If (!Exists(R_PML)) R_PML=0.4; EndIf
If (!Exists(R_TR)) R_TR=0.1; EndIf
        // Ordinate (theta coordinate in rad)
If (!Exists(phi)) phi =90*Pi/180; EndIf // corner angle
        // Offset
If (!Exists(x_ofst)) x_ofst = 0; EndIf
If (!Exists(y_ofst)) y_ofst = 0; EndIf
        // Mesh
If (!Exists(Nz)) Nz=10; EndIf // No. of elements in z-direction
If (!Exists(Ntheta)) Ntheta=20; EndIf // No. of elements in theta-direction    
		// Macro directory (no trailing '/')
If (!Exists(MACRO_DIR)) MACRO_DIR="."; EndIf

    // ---- Corner geometry
name = "corner";
    // Abscissa
x_input[] = {Log(R_TR), Log(R_PML), Log(R)};
    // Angle
y_input[] = {-Pi, -phi/2, phi/2, Pi};

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

    // Mesh in theta-direction
For i In {0:(#y_input[]-2)}
    N_mesh_y[i] = Ceil(Ntheta*(y_input[i+1]-y_input[i])/(y_input[#y_input[]-1]-y_input[0]))+1;
EndFor

    // Definition of physical entities
    // Set to '0' to disable default physical entities
Bool_define_line_physical_entities = 1;
Bool_define_surface_physical_entities = 0;

//MACRO_DIR = "/home/bino/shared/Plasmonics/mesh";
//MACRO_DIR = ".";
//Include StrCat(MACRO_DIR,"/Macro-Rectangle-Array.geo");
Include StrCat(CurrentDirectory,"Macro-Rectangle-Array.geo");
//Include "Macro-Rectangle-Array.geo";


    // Physical entities
Physical Surface("omega-m") = {s~{4}};
Physical Surface("omega-d") = {s~{1},s~{7}};
Physical Surface("omega-m-pml") = {s~{3}};
Physical Surface("omega-d-pml") = {s~{0},s~{6}};
