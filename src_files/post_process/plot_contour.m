function[] = plot_contour(phi)

global Elements;global Boundaries;

if length(phi)==length(Elements.centroid)
    pos = Elements.centroid;
else
    pos = [Elements.centroid;Boundaries.loc_bound];
end

Nodes = Elements.Nodes.coord;
interpol = scatteredInterpolant(pos,phi);

phi_nodal = interpol(Nodes);

patch('Vertices',Nodes,'Faces',Elements.faces.nodes,'EdgeColor','None','FaceVertexCData',phi_nodal,'FaceColor','interp');
colorbar

end