%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the convection term of the Navier-stokes over the
% entire domain, i.e. integral(del(uu) dv)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[conv] = convection(Boundaries,Elements,u,v,p)

a = Elements.faces;
[~,nd] = size(a.nodes);
centroids = Elements.centroid;
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;
[ugrad,vgrad,~] = gradient_loop(Elements,Boundaries,u,v,p);


j=1;
x_conv = zeros(length(centroids),1);
y_conv = x_conv;
for i=1:nd
    
    nx = a.normal(:,j);
    ny = a.normal(:,j+1);
    ds = a.area(:,i);
    faces = a.mid(:,[j,j+1]);
    neighb_pos = a.neighb(:,i);
    neighb_bound = a.neighb_bound(:,i);
    neighb = centroids(neighb_pos,:);
%     fx = vecnorm(faces-neighb,2,2)./vecnorm(centroids-neighb,2,2);
%     fx(abs(fx)==Inf)=0;
    bd = a.bound_flag(:,i);

    r = faces - 0.5*(centroids+neighb);
    uf = 0.5*(u+u(neighb_pos)) + 0.5*dot(ugrad+ugrad(neighb_pos),r,2);
    vf = 0.5*(v+v(neighb_pos)) + 0.5*dot(vgrad+vgrad(neighb_pos),r,2);
    
    uf(bd==1) = u_bound(neighb_bound(bd==1));
    vf(bd==1) = v_bound(neighb_bound(bd==1));
    
    x_conv = x_conv + uf.*(uf.*nx + vf.*ny).*ds;
    y_conv = y_conv + vf.*(uf.*nx + vf.*ny).*ds;
    
    
    j=j+2;

end

conv = [x_conv,y_conv];

end
