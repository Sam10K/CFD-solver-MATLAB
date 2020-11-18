function[cont,dt] = continuity(Elements,varargin)

a = Elements.faces;
[~,nd] = size(a.nodes);

j=1;
cont = zeros(length(Elements.centroid),1);
abs_cont = zeros(length(Elements.centroid),1);

for i=1:nd
    
    nx = a.normal(:,j);
    ny = a.normal(:,j+1);
    ds = a.area(:,i);
    uf = Elements.faces.uf(:,i);
    vf = Elements.faces.vf(:,i);
    
    cont = cont + (uf.*nx + vf.*ny).*ds;
    abs_cont = abs_cont + abs((uf.*nx + vf.*ny).*ds);
    
    j=j+2;

end

if ~isempty(varargin)
    CFL = varargin{1};
    dt = 2*CFL*Elements.volume./abs_cont;
end

end
