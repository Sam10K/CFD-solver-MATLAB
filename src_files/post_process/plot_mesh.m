function[] = plot_mesh()

global Elements;

Nodes = Elements.Nodes.coord;

patch('Vertices',Nodes,'Faces',Elements.faces.nodes,'EdgeColor','k','FaceColor','None');


end