function[Elements] = face_interp(Elements,Boundaries,u,v)

a = Elements.faces;
[~,nd] = size(a.nodes);
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;

centroids = [Elements.centroid];

j=1;

for i=1:nd

    nx = a.normal(:,j);
    ny = a.normal(:,j+1);
    ds = a.area(:,i);
    faces = a.mid(:,[j,j+1]);
    S = [ds.*nx,ds.*ny];
    neighb_pos = a.neighb(:,i);
    neighb = centroids(neighb_pos,:);
    bd = logical(a.bound_flag(:,i));    
    
    neighb_bound = a.neighb_bound(:,i);        
    
%     close
%     quiver(faces(bd,1),faces(bd,2),nx(bd),ny(bd))
%     patch('Vertices',[Elements.Nodes.coord(:,1),Elements.Nodes.coord(:,2)],'Faces',[Elements.faces.nodes],'FaceColor','None')
%     pause

    fx = vecnorm(faces-neighb,2,2)./vecnorm(centroids-neighb,2,2);
    fx(abs(fx)==Inf)=0;fx(isnan(fx)) = 0;
    Elements.faces.uf(:,i) = fx.*u + (1-fx).*u(neighb_pos);
    Elements.faces.vf(:,i) = fx.*v + (1-fx).*v(neighb_pos);
    Elements.faces.uf(bd,i) = u_bound(neighb_bound(bd));
    Elements.faces.vf(bd,i) = v_bound(neighb_bound(bd));
    
    j=j+2;
    
end



end