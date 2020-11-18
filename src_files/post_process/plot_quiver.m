function[] = plot_quiver(u,v)

global Elements;

quiver(Elements.centroid(:,1),Elements.centroid(:,2),u,v,'Color','k');


end