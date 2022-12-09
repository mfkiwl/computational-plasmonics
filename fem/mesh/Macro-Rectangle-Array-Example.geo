/*
* This geometry file demonstrate the addition of a corner region meshed in Euler coordinates to a main geometry.
* The steps are:
*   1) The main geometry is constructed, truncated by a disk of radius R.
*   2) The corner geometry is added using an 'Import' statement.
*   3) Periodicity conditions are defined, to map the mesh at the truncation radius R to the corner mesh.
*/
SetFactory("OpenCASCADE");

    // ---- Input parameters

R_O = 1; // outer radius
R = 0.5*R_O; // corner radius
R_PML = 0.2*R_O; // PML radius
R_TR = 0.01*R_O; // Truncation radius
    
xc = 0; yc = 0; // corner coordinates

theta_i = Pi/3; // inner angle
phi = Pi/2; // corner angle
theta_o = 2*phi-theta_i; // outer angle
        // mesh
Mesh.CharacteristicLengthMin = 0.05;
Mesh.CharacteristicLengthMax = 0.05;
//Mesh.ElementOrder = 2;

    // ---- Geometry
        // outer circle
l_circ_out = newl; Circle(l_circ_out) = {xc, yc, 0, R_O, 0, 2*Pi};
        // inner circle (top-half)
l_circ_1 = newl; Circle(l_circ_1) = {xc, yc, 0, R, -theta_i/2, theta_i/2};
l_circ_2 = newl; Circle(l_circ_2) = {xc, yc, 0, R, theta_i/2, phi/2};
l_circ_3 = newl; Circle(l_circ_3) = {xc, yc, 0, R, phi/2, theta_o/2};
l_circ_4 = newl; Circle(l_circ_4) = {xc, yc, 0, R, theta_o/2, Pi};
        // inner circle (bottom-half)
l_circ_5 = newl; Circle(l_circ_5) = {xc, yc, 0, R, Pi, -theta_o/2};
l_circ_6 = newl; Circle(l_circ_6) = {xc, yc, 0, R, -theta_o/2, -phi/2};
l_circ_7 = newl; Circle(l_circ_7) = {xc, yc, 0, R, -phi/2, -theta_i/2};

l_circ_out_loop = newl; Curve Loop(l_circ_out_loop) = {l_circ_out};
l_circ_loop = newl; Curve Loop(l_circ_loop) = {l_circ_1,l_circ_2,l_circ_3,l_circ_4,l_circ_5,l_circ_6,l_circ_7};

s_disk = news; Plane Surface(s_disk) = {l_circ_out_loop,l_circ_loop};
    // ---- Physical entities
Physical Curve ("Disk-Outer-Boundary") = {l_circ_out};
Physical Curve ("Disk-Inner-Boundary-1") = {l_circ_1};
Physical Curve ("Disk-Inner-Boundary-2") = {l_circ_2};
Physical Curve ("Disk-Inner-Boundary-3") = {l_circ_3};
Physical Curve ("Disk-Inner-Boundary-4") = {l_circ_4};
Physical Curve ("Disk-Inner-Boundary-5") = {l_circ_5};
Physical Curve ("Disk-Inner-Boundary-6") = {l_circ_6};
Physical Curve ("Disk-Inner-Boundary-7") = {l_circ_7};
Physical Surface ("Disk") = {s_disk};

    // ---- Corner geometry
name = "corner";
    // Abscissa
x_input[] = { Log(R_TR), Log(R_PML), Log(R)};
    // Angle
y_input[] = { -Pi, -theta_o/2, -phi/2,-theta_i/2,theta_i/2,phi/2,theta_o/2,Pi};
        // Offset
x_ofst = -R_O;
y_ofst = -R_O-Pi;

For i In {0:(#x_input[]-1)}
    x_input[i] = x_input[i] + x_ofst;
EndFor

For i In {0:(#y_input[]-1)}
    y_input[i] = y_input[i] + y_ofst;
EndFor
    // Mesh in z-direction
Nz = 20; // total number of element in z direction
For i In {0:(#x_input[]-2)}
    N_mesh_x[i] = Ceil(Nz*(x_input[i+1]-x_input[i])/(x_input[#x_input[]-1]-x_input[0]))+1;
EndFor

    // Definition of physical entities
    // Set to '0' to disable default physical entities
Bool_define_line_physical_entities = 1;
Bool_define_surface_physical_entities = 1;

Include "Macro-Rectangle-Array.geo";

    // Periodicity disk boundary -> Right boundary
l_src[] = {l_circ_5, l_circ_6, l_circ_7, l_circ_1, l_circ_2, l_circ_3, l_circ_4};
Call ConnectSrcToRightBnd;
