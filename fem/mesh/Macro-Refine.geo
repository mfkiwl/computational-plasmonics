// Macro refine
    // Refine by splitting (reduces order)
order = Mesh.ElementOrder;
RefineMesh;
    // Restore original mesh order
SetOrder order;
