function[u,v,p,Boundaries,Elements] = initialize(BC,Elements,Boundaries)

u = BC.internal.U.value; if size(u,1)==1 u=u*ones(length(Elements.centroid),1);end
v = BC.internal.V.value; if size(v,1)==1 v=v*ones(length(Elements.centroid),1);end
p = BC.internal.P.value; if size(p,1)==1 p=p*ones(length(Elements.centroid),1);end

[Boundaries] = boundary_conditions(Boundaries,u,v,p,BC);
centroids = [Elements.centroid];
loc_bound = [Boundaries.loc_bound];
ubound_type=Boundaries.ubound_type;
pbound_type=Boundaries.pbound_type;
a = Elements.faces;
[~,nd] = size(a.nodes);
Elements.faces.neighb_bound = zeros(length(centroids),nd);
Elements.faces.ubound_type = zeros(length(centroids),nd);
Elements.faces.pbound_type = zeros(length(centroids),nd);

for i=1:length(centroids)
    j=1;
    for k=1:nd
    
       face = a.mid(i,[j,j+1]);
       dist = (face(1)-loc_bound(:,1)).^2 + (face(2)-loc_bound(:,2)).^2;
       [d,pos] = min(dist);
       if d==0
           Elements.faces.neighb_bound(i,k) = pos;
           Elements.faces.ubound_type(i,k) = ubound_type(pos);
           Elements.faces.pbound_type(i,k) = pbound_type(pos);
       else
           Elements.faces.neighb_bound(i,k) = 0;
       end
       j=j+2;
    end
   
end

end