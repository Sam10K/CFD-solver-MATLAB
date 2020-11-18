function[ugrad,vgrad,pgrad] = gradient_loop_lsq(Elements,Boundaries,u,v,p)

a = Elements.faces;
[~,nd] = size(a.nodes);
centroids = Elements.centroid;
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;
p_bound = Boundaries.p_bound;

j=1;

a11=zeros(size(u));a22=a11;a12=a11;a21=a11;
b1_u=a11;b2_u=a11;
b1_v=a11;b2_v=a11;
b1_p=a11;b2_p=a11;
dummy=[];
for i=1:nd
    
    faces = a.mid(:,[j,j+1]);
    neighb_pos = a.neighb(:,i);
    neighb_bound = a.neighb_bound(:,i);
    neighb = centroids(neighb_pos,:);
    un = u(neighb_pos);
    vn = v(neighb_pos);
    pn = p(neighb_pos);
    
    bd = logical(a.bound_flag(:,i));
    neighb(bd,:) = faces(bd,:);
    un(bd) = u_bound(neighb_bound(bd));
    vn(bd) = v_bound(neighb_bound(bd));
    pn(bd) = p_bound(neighb_bound(bd));
    
    delta_x = neighb(:,1)-centroids(:,1);
    delta_y = neighb(:,2)-centroids(:,2);
    w = 1./vecnorm(neighb-centroids,2,2);
    
    delta_u = un - u;
    delta_v = vn - v;
    delta_p = pn - p;
    
    
    a11 = a11 + w.*delta_x.*delta_x;
    a12 = a12 + w.*delta_x.*delta_y;
    a21 = a21 + w.*delta_y.*delta_x;
    a22 = a22 + w.*delta_y.*delta_y;
    
    b1_u = b1_u + w.*delta_x.*delta_u;
    b2_u = b2_u + w.*delta_y.*delta_u;
    
    b1_v = b1_v + w.*delta_x.*delta_v;
    b2_v = b2_v + w.*delta_y.*delta_v;
    
    b1_p = b1_p + w.*delta_x.*delta_p;
    b2_p = b2_p + w.*delta_y.*delta_p;
    
    
    j=j+2;
    
   

end

dudy = (b2_u.*a11 - a21.*b1_u) ./ (a22.*a11 - a21.*a12);
dudx = (b1_u - a12.*dudy)./a11;

dvdy = (b2_v.*a11 - a21.*b1_v) ./ (a22.*a11 - a21.*a12);
dvdx = (b1_v - a12.*dvdy)./a11;

dpdy = (b2_p.*a11 - a21.*b1_p) ./ (a22.*a11 - a21.*a12);
dpdx = (b1_p - a12.*dpdy)./a11;

ugrad = [dudx,dudy];
vgrad = [dvdx,dvdy];
pgrad = [dpdx,dpdy];



end