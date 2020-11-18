function[A_u,B_u,A_v,B_v] = H_Ap_assemble(Boundaries,Elements,Properties,u,v,p,Schemes,varargin)

nu = Properties.nu;
rho = Properties.rho;

a = Elements.faces;
[~,nd] = size(a.nodes);
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;

centroids = [Elements.centroid];
[ugrad,vgrad,~] = gradient_loop(Elements,Boundaries,u,v,p,Schemes);

j=1;

A_u = sparse(size(u,1),size(u,1)); B_u = zeros(size(u));
A_v = sparse(size(u,1),size(u,1)); B_v = zeros(size(u)); 

for i=1:nd

%% Diffusion 
    nx = a.normal(:,j);
    ny = a.normal(:,j+1);
    ds = a.area(:,i);
    faces = a.mid(:,[j,j+1]);
    S = [ds.*nx,ds.*ny];
    neighb_pos = a.neighb(:,i);
    neighb = centroids(neighb_pos,:);
    d = (neighb-centroids);
    delta = d./dot(d,S,2) .* vecnorm(S,2,2).^2;
    bd = logical(a.bound_flag(:,i));
    ubound_type = a.ubound_type(:,i);
    bdzg = ubound_type==2;
    bdw = ubound_type==3;
    neighb_bound = a.neighb_bound(:,i);
    un = u(neighb_pos);
    vn = v(neighb_pos);
    db = (faces - centroids);
    dn = dot(db,S,2).*S./(S(:,1).^2 + S(:,2).^2);
    un(bd,:)=  u_bound(neighb_bound(bd));
    vn(bd,:)=  v_bound(neighb_bound(bd));
    d(bd,:)=dn(bd,:);
    delta(bd,:)=S(bd,:);
    k = S - delta;

    del_u = grad_interp(faces,centroids,neighb,neighb_pos,u,ugrad,bd,S,neighb_bound,u_bound);
    del_v = grad_interp(faces,centroids,neighb,neighb_pos,v,vgrad,bd,S,neighb_bound,v_bound);
    del_u(bdzg,:) = 0;
    del_v(bdzg,:) = 0;
    
    Del_u_velocity = [del_u(:,1),del_v(:,1)];% First row of Del_u tensor
    Del_v_velocity = [del_u(:,2),del_v(:,2)];% Second row of Del_u tensor
    
    uf = Elements.faces.uf(:,i);
    vf = Elements.faces.vf(:,i);
    F = (uf.*nx + vf.*ny).*ds;
    [~,~,b_u] = interpol_coeff(faces,centroids,neighb,neighb_pos,u,ugrad,bd,S,neighb_bound,u_bound,F,Schemes);
    [~,~,b_v] = interpol_coeff(faces,centroids,neighb,neighb_pos,v,vgrad,bd,S,neighb_bound,v_bound,F,Schemes);
    
    % Eq form : ac*uc + af*uf = b
    ap_xdiff = nu*vecnorm(delta,2,2)./vecnorm(d,2,2);
    an_xdiff = -nu*vecnorm(delta,2,2)./vecnorm(d,2,2);
    b_xdiff = (nu*dot(Del_u_velocity,k,2) + nu*dot(del_u,S,2));
    ap_ydiff = ap_xdiff;
    an_ydiff = an_xdiff;
    b_ydiff = (nu*dot(Del_v_velocity,k,2) + nu*dot(del_v,S,2));
    
    ap_xconv = max(F,0);
    an_xconv = -max(-F,0);
    b_xconv = -F.*b_u;
    ap_yconv = ap_xconv;
    an_yconv = an_xconv;
    b_yconv = -F.*b_v;
    
    % Boundary face contributions
    
    %% Wall boundary
    
    ap_xdiff_dummy = nu*(1-nx.^2).*vecnorm(S,2,2)./vecnorm(d,2,2);
    ap_xdiff(bdw) = ap_xdiff_dummy(bdw);
    an_xdiff(bdw) = 0;
    b_bdw_dummy = nu*vecnorm(S,2,2)./(vecnorm(d,2,2)) .* ...
                  (un.*(1-nx.^2)+(v-vn).*nx.*ny);             
    b_xdiff(bdw) = b_bdw_dummy(bdw);
    
    ap_ydiff_dummy = nu*(1-ny.^2).*vecnorm(S,2,2)./vecnorm(d,2,2);
    ap_ydiff(bdw) = ap_ydiff_dummy(bdw);
    an_ydiff(bdw) = 0;
    b_bdw_dummy = nu*vecnorm(S,2,2)./(vecnorm(d,2,2)) .* ...
                  (vn.*(1-ny.^2)+(u-un).*nx.*ny);
    b_ydiff(bdw) = b_bdw_dummy(bdw);
    
    ap_xconv(bdw) = 0;
    an_xconv(bdw) = 0;
    b_xconv(bdw) = 0;
    
    ap_yconv(bdw) = 0;
    an_yconv(bdw) = 0;
    b_yconv(bdw) = 0;              
    
    %% Zero gradient boundary
    
    ap_xdiff(bdzg) = 0;
    an_xdiff(bdzg) = 0;
    b_xdiff(bdzg) = 0;
    
    ap_ydiff(bdzg) = 0;    
    an_ydiff(bdzg) = 0;    
    b_ydiff(bdzg) = 0;
    
    ap_xconv(bdzg) = F(bdzg);
    an_xconv(bdzg) = 0;
    b_xconv(bdzg)=0;
    
    ap_yconv(bdzg) = F(bdzg);
    an_yconv(bdzg) = 0;
    b_yconv(bdzg)=0;
    
    %% Specified velocity boundary
    bd_other = logical(bd.* ~bdzg.* ~bdw); % Boundaries other than zg and wall
    ap_xdiff(bd_other) = ap_xdiff(bd_other);    
    b_xdiff(bd_other) = -an_xdiff(bd_other).*un(bd_other) + b_xdiff(bd_other);
    an_xdiff(bd_other) = 0;
    
    ap_ydiff(bd_other) = ap_ydiff(bd_other);    
    b_ydiff(bd_other) = -an_ydiff(bd_other).*vn(bd_other) + b_ydiff(bd_other);
    an_ydiff(bd_other) = 0;
    
    ap_xconv(bd_other) = 0;
    an_xconv(bd_other) = 0;
    b_xconv(bd_other) = -un(bd_other).*F(bd_other);
    
    ap_yconv(bd_other) = 0;
    an_yconv(bd_other) = 0;
    b_yconv(bd_other) = -vn(bd_other).*F(bd_other);
   
    %% Final sparse matrices assembly
        
    a_c_u = ap_xconv + ap_xdiff; a_f_u = an_xconv + an_xdiff; b_u = b_xconv + b_xdiff;
    a_c_v = ap_yconv + ap_ydiff; a_f_v = an_yconv + an_ydiff; b_v = b_yconv + b_ydiff;
        
    A_u = A_u + rho*sparse([1:size(u,1)]',[1:size(u,1)]',a_c_u,size(u,1),size(u,1)) + ...
                rho*sparse([1:size(u,1)]',neighb_pos,a_f_u,size(u,1),size(u,1));
    
    A_v = A_v + rho*sparse([1:size(u,1)]',[1:size(u,1)]',a_c_v,size(u,1),size(u,1)) + ...
                rho*sparse([1:size(u,1)]',neighb_pos,a_f_v,size(u,1),size(u,1));
    
    B_u = B_u + rho*b_u;
    B_v = B_v + rho*b_v;
                
    j=j+2;
     
end

if ~isempty(varargin)
    
    a=1;b=0;
    
    dt=varargin{1};u_old=varargin{2};v_old=varargin{3};
    conv_n=varargin{4};diff_n=varargin{5};
    
    diag_col = rho*Elements.volume/dt;
    A_u = a*A_u + spdiags(diag_col,0,size(u,1),size(u,1));
    A_v = a*A_v + spdiags(diag_col,0,size(u,1),size(u,1));

    B_u = rho*u_old.*Elements.volume/dt + a*B_u + b*rho*(diff_n(:,1)-conv_n(:,1));
    B_v = rho*v_old.*Elements.volume/dt + a*B_v + b*rho*(diff_n(:,2)-conv_n(:,2));
    
end

end
