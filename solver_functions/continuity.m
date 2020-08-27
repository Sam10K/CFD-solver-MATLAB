%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate continuity over the entire domain and the time step
% size required for the next time step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[cont,dt] = continuity(Boundaries,Elements,u,v,CFL)

a = Elements.faces;
[~,nd] = size(a.nodes);
centroids = Elements.centroid;
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;

j=1;
cont = zeros(length(centroids),1);
abs_cont=cont;

for i=1:nd
    
    nx = a.normal(:,j);
    ny = a.normal(:,j+1);
    ds = a.area(:,i);
    faces = a.mid(:,[j,j+1]);
    neighb_pos = a.neighb(:,i);
    neighb_bound = a.neighb_bound(:,i);
    neighb = centroids(neighb_pos,:);
    fx = vecnorm(faces-neighb,2,2)./vecnorm(centroids-neighb,2,2);
    fx(abs(fx)==Inf)=0;
    uf = fx.*u + (1-fx).*u(neighb_pos);
    vf = fx.*v + (1-fx).*v(neighb_pos);
    bd = a.bound_flag(:,i);
    uf(bd==1) = u_bound(neighb_bound(bd==1));
    vf(bd==1) = v_bound(neighb_bound(bd==1));
    
    cont = cont + (uf.*nx + vf.*ny).*ds;
    
    abs_cont = abs_cont + abs((uf.*nx + vf.*ny).*ds);
    
    j=j+2;

end

dt = 2*CFL*Elements.volume./abs_cont;


end
