function[Elements] = rhie_chow(Elements,Boundaries,u,v,p,Ap,alpha_u,Schemes)

a = Elements.faces;
[~,nd] = size(a.nodes);
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;
p_bound = Boundaries.p_bound;

centroids = [Elements.centroid];
[ugrad,vgrad,pgrad] = gradient_loop(Elements,Boundaries,u,v,p,Schemes);

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
    
    uf_dummy = Elements.faces.uf(:,i);
    vf_dummy = Elements.faces.vf(:,i);
%     F = (uf_dummy.*nx + vf_dummy.*ny).*ds;
% 
%     [ap_u,an_u,b_u] = interpol_coeff(faces,centroids,neighb,neighb_pos,u,ugrad,bd,S,neighb_bound,u_bound,F);
%     uf = ap_u.*u + an_u.*u(neighb_pos) + b_u;
% 
%     [ap_v,an_v,b_v] = interpol_coeff(faces,centroids,neighb,neighb_pos,v,vgrad,bd,S,neighb_bound,v_bound,F);
%     vf = ap_v.*v + an_v.*v(neighb_pos) + b_v;
    
    fx = vecnorm(faces-neighb,2,2)./vecnorm(centroids-neighb,2,2);
    fx(abs(fx)==Inf)=0;
    fx(isnan(fx)) = 0;
%     
    uf = fx.*u + (1-fx).*u(neighb_pos);
    vf = fx.*v + (1-fx).*v(neighb_pos);
    uf(bd) = u_bound(neighb_bound(bd));
    vf(bd) = v_bound(neighb_bound(bd));
    
    del_p_bar = [fx.*pgrad(:,1) + (1-fx).*pgrad(neighb_pos,1),...
                 fx.*pgrad(:,2) + (1-fx).*pgrad(neighb_pos,2)];
    del_p_bar(bd,:) = pgrad(bd,:);
    
    del_p = grad_interp(faces,centroids,neighb,neighb_pos,p,pgrad,bd,S,neighb_bound,p_bound);
    
    Ap_f = [fx.*Ap(:,1) + (1-fx).*Ap(neighb_pos,1),...
            fx.*Ap(:,2) + (1-fx).*Ap(neighb_pos,2)];
    Ap_f(bd,:) = Ap(bd,:);
    
    uf = uf - 1./Ap_f(:,1) .* (del_p(:,1)-del_p_bar(:,1)) + (1-alpha_u)*(uf_dummy-uf);
    vf = vf - 1./Ap_f(:,2) .* (del_p(:,2)-del_p_bar(:,2)) + (1-alpha_u)*(vf_dummy-vf);
    
%     uf(bd) = u_bound(neighb_bound(bd));
%     vf(bd) = v_bound(neighb_bound(bd));
    
    Elements.faces.uf(:,i) = uf;
    Elements.faces.vf(:,i) = vf;
    
    j=j+2;
    
end



end