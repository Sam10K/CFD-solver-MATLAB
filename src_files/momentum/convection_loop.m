function[conv] = convection_loop(Boundaries,Elements,u,v,p,Schemes)

a = Elements.faces;
[~,nd] = size(a.nodes);
centroids = Elements.centroid;
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;
[ugrad,vgrad,~] = gradient_loop(Elements,Boundaries,u,v,p,Schemes);


j=1;
x_conv = zeros(length(centroids),1);
y_conv = x_conv;
for i=1:nd
    
    nx = a.normal(:,j);
    ny = a.normal(:,j+1);
    ds = a.area(:,i);
    S = [ds.*nx,ds.*ny];
    faces = a.mid(:,[j,j+1]);
    neighb_pos = a.neighb(:,i);
    neighb_bound = a.neighb_bound(:,i);
    neighb = centroids(neighb_pos,:);

    bd = a.bound_flag(:,i);
    
    uf_dummy = Elements.faces.uf(:,i);
    vf_dummy = Elements.faces.vf(:,i);
    
    F = (uf_dummy.*nx + vf_dummy.*ny).*ds;
    
    [ap,an,b] = interpol_coeff(faces,centroids,neighb,neighb_pos,u,ugrad,bd,S,neighb_bound,u_bound,F,Schemes);
    uf = ap.*u + an.*u(neighb_pos) + b;
    
    [ap,an,b] = interpol_coeff(faces,centroids,neighb,neighb_pos,v,vgrad,bd,S,neighb_bound,v_bound,F,Schemes);
    vf = ap.*v + an.*v(neighb_pos) + b;
    
    uf(bd==1) = u_bound(neighb_bound(bd==1));
    vf(bd==1) = v_bound(neighb_bound(bd==1));
    
    x_conv = x_conv + uf.*(uf.*nx + vf.*ny).*ds;
    y_conv = y_conv + vf.*(uf.*nx + vf.*ny).*ds;
    
    
    j=j+2;

end

conv = [x_conv,y_conv];

end
