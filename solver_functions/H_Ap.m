%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate H and Ap vectors for solving the momentum equations
% and find the velocities ustar, vstar from the equation :
% Ap*U = H - del(p) ==> U = H./Ap - del(p)./Ap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function[H,Ap] = H_Ap(Boundaries,Elements,Properties,u,v,p,conv_n,diff_n,dt,u_old,v_old)

nu = Properties.nu;
rho = Properties.rho;

a = Elements.faces;
[~,nd] = size(a.nodes);
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;

centroids = [Elements.centroid];
[Ugrad,Vgrad,~] = gradient(Elements,Boundaries,u,v,p);
dudx_p = Ugrad(:,1);
dudy_p = Ugrad(:,2);
dvdx_p = Vgrad(:,1);
dvdy_p = Vgrad(:,2);

j=1;

Ap_xdiff = zeros(size(u)); B_xdiff = zeros(size(u)); BC_xdiff = zeros(size(u));
Ap_ydiff = zeros(size(u)); B_ydiff = zeros(size(u)); BC_ydiff = zeros(size(u));
Ap_xconv = zeros(size(u)); B_xconv = zeros(size(u)); BC_xconv = zeros(size(u));
Ap_yconv = zeros(size(u)); B_yconv = zeros(size(u)); BC_yconv = zeros(size(u));

for i=1:nd

%% Diffusion 
    nx = a.normal(:,j);
    ny = a.normal(:,j+1);
    ds = a.area(:,i);
    faces = a.mid(:,[j,j+1]);
    S = [ds.*nx,ds.*ny];
    neighb_pos = a.neighb(:,i);
    neighb = centroids(neighb_pos,:);
    d = neighb-centroids;
    delta = d./dot(d,S,2) .* vecnorm(S,2,2).^2;%d.*vecnorm(S,2,2)./vecnorm(d,2,2);
    k = S - delta;
    un = u(neighb_pos);
    vn = v(neighb_pos);
    up = u;
    vp = v;
    fx = vecnorm(faces-neighb,2,2)./vecnorm(centroids-neighb,2,2);
    fx(abs(fx)==Inf)=0;
    del_u = [dudx_p.*fx + (1-fx).*dudx_p(neighb_pos),dudy_p.*fx + (1-fx).*dudy_p(neighb_pos)];
    del_v = [dvdx_p.*fx + (1-fx).*dvdx_p(neighb_pos),dvdy_p.*fx + (1-fx).*dvdy_p(neighb_pos)];
    bd = Elements.faces.bound_flag(:,i);
    db = faces - centroids;
    dn = dot(db,S,2).*S./(S(:,1).^2 + S(:,2).^2);
    k(bd==1,:)=0;
    neighb_bound = a.neighb_bound(:,i);
    un(bd==1,:)=  u_bound(neighb_bound(bd==1));
    vn(bd==1,:)=  v_bound(neighb_bound(bd==1));
    d(bd==1,:)=dn(bd==1,:);
    delta(bd==1,:)=S(bd==1,:);
    e = d./vecnorm(d,2,2);

    % Refer page 289 of uFVM book
    del_u = [del_u(:,1) + ((un-up)./vecnorm(d,2,2) - dot(del_u,e,2)).*e(:,1),...
              del_u(:,2) + ((un-up)./vecnorm(d,2,2) - dot(del_u,e,2)).*e(:,2)];
    del_v = [del_v(:,1) + ((vn-vp)./vecnorm(d,2,2) - dot(del_v,e,2)).*e(:,1),...
              del_v(:,2) + ((vn-vp)./vecnorm(d,2,2) - dot(del_v,e,2)).*e(:,2)];

    % Matrices for xdiff
    ap_xdiff = -nu*vecnorm(delta,2,2)./vecnorm(d,2,2);
    b_xdiff = nu*(vecnorm(delta,2,2).*un./vecnorm(d,2,2) + dot(k,del_u,2));
    bc_xdiff = zeros(size(u));
    xdiff_dummy = nu*(vecnorm(delta,2,2).*(un-up)./vecnorm(d,2,2) + dot(k,del_u,2));
    bc_xdiff(bd==1) = -(ap_xdiff(bd==1).*up(bd==1)+b_xdiff(bd==1)) + xdiff_dummy(bd==1);

    % Matrices for ydiff
    ap_ydiff = -nu*vecnorm(delta,2,2)./vecnorm(d,2,2);
    b_ydiff = nu*(vecnorm(delta,2,2).*vn./vecnorm(d,2,2) + dot(k,del_v,2));
    bc_ydiff = zeros(size(u));
    ydiff_dummy = nu*(vecnorm(delta,2,2).*(vn-vp)./vecnorm(d,2,2) + dot(k,del_v,2));
    bc_ydiff(bd==1) = -(ap_ydiff(bd==1).*vp(bd==1)+b_ydiff(bd==1)) + ydiff_dummy(bd==1);


    Ap_xdiff = Ap_xdiff + ap_xdiff;
    B_xdiff = B_xdiff + b_xdiff;
    BC_xdiff = BC_xdiff + bc_xdiff;

    Ap_ydiff = Ap_ydiff + ap_ydiff;
    B_ydiff = B_ydiff + b_ydiff;
    BC_ydiff = BC_ydiff + bc_ydiff;

%% Convection

    r = faces - 0.5*(centroids+neighb);
    uf = 0.5*(u+u(neighb_pos)) + 0.5*dot(Ugrad+Ugrad(neighb_pos),r,2);
    vf = 0.5*(v+v(neighb_pos)) + 0.5*dot(Vgrad+Vgrad(neighb_pos),r,2);
    
    uf(bd==1) = u_bound(neighb_bound(bd==1));
    vf(bd==1) = v_bound(neighb_bound(bd==1));
    
    %Matrices for xconv
    ap_xconv = 0.5*(uf.*nx + vf.*ny).*ds;
    b_xconv = 0.5*u(neighb_pos).*(uf.*nx + vf.*ny).*ds + ...
        0.5*dot(Ugrad+Ugrad(neighb_pos),r,2).*(uf.*nx + vf.*ny).*ds;
    
    xconv_dummy = uf.*(uf.*nx + vf.*ny).*ds;
    
    bc_xconv = zeros(size(u));
    
    bc_xconv(bd==1) = -(ap_xconv(bd==1).*u(bd==1) + b_xconv(bd==1))...
                        + xconv_dummy(bd==1);
    
    %Matrices for yconv                
    ap_yconv = 0.5*(uf.*nx + vf.*ny).*ds;
    b_yconv = 0.5*v(neighb_pos).*(uf.*nx + vf.*ny).*ds + ...
        0.5*dot(Vgrad+Vgrad(neighb_pos),r,2).*(uf.*nx + vf.*ny).*ds;
    
    yconv_dummy = vf.*(uf.*nx + vf.*ny).*ds;
    
    bc_yconv = zeros(size(u));
    
    bc_yconv(bd==1) = -(ap_yconv(bd==1).*v(bd==1) + b_yconv(bd==1))...
                        + yconv_dummy(bd==1);
    

    
    Ap_xconv = Ap_xconv + ap_xconv;
    B_xconv = B_xconv + b_xconv;
    BC_xconv = BC_xconv + bc_xconv;
    
    Ap_yconv = Ap_yconv + ap_yconv;
    B_yconv = B_yconv + b_yconv;
    BC_yconv = BC_yconv + bc_yconv;


j=j+2;
     
end

x_conv_n = conv_n(:,1);y_conv_n = conv_n(:,2);
x_diff_n = diff_n(:,1);y_diff_n = diff_n(:,2);

Ap_x = (1+0.5*dt*Ap_xconv./Elements.volume - 0.5*dt*Ap_xdiff./Elements.volume)*rho/dt;
      
H_x = u_old*rho/dt + 0.5*(-B_xconv-BC_xconv-x_conv_n + B_xdiff+BC_xdiff+x_diff_n)*rho./Elements.volume;
       
Ap_y = (1+0.5*dt*Ap_yconv./Elements.volume - 0.5*dt*Ap_ydiff./Elements.volume)*rho/dt;
       
H_y = v_old*rho/dt + 0.5*(-B_yconv-BC_yconv-y_conv_n + B_ydiff+BC_ydiff+y_diff_n)*rho./Elements.volume;

H = [H_x,H_y];
Ap = [Ap_x,Ap_y];

end
