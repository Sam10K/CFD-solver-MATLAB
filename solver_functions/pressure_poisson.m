%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to solve the pressure poisson equation to obtain the pressure across 
% the domain. Includes Rhie-Chow interpolation to the face velocities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[p_centre] = pressure_poisson(Boundaries,Elements,Properties,u,v,p,dt)

rho = Properties.rho;

a = Elements.faces;
[~,nd] = size(a.nodes);
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;
p_bound = Boundaries.p_bound;

centroids = [Elements.centroid];
[Ugrad,Vgrad,Pgrad] = gradient_loop(Elements,Boundaries,u,v,p);
dpdx_p = Pgrad(:,1);
dpdy_p = Pgrad(:,2);

j=1;
rhs = zeros(length(centroids),1);
k_delp = rhs;pn_delta_d = rhs; delta_d = rhs;

for i=1:nd
    

nx = a.normal(:,j);
ny = a.normal(:,j+1);
ds = a.area(:,i);
faces = a.mid(:,[j,j+1]);
S = [ds.*nx,ds.*ny];
neighb_pos = a.neighb(:,i);
neighb = centroids(neighb_pos,:);
fx = vecnorm(faces-neighb,2,2)./vecnorm(centroids-neighb,2,2);
fx(abs(fx)==Inf)=0;
% uf = fx.*u + (1-fx).*u(neighb_pos);
% vf = fx.*v + (1-fx).*v(neighb_pos);
r = faces-0.5*(centroids+neighb);
uf = 0.5*(u+u(neighb_pos)) + 0.5*(dot(Ugrad+Ugrad(neighb_pos),r,2));
vf = 0.5*(v+v(neighb_pos)) + 0.5*(dot(Vgrad+Vgrad(neighb_pos),r,2));
d = neighb-centroids;
delta = d./dot(d,S,2) .* vecnorm(S,2,2).^2;%delta = d.*vecnorm(S,2,2)./vecnorm(d,2,2);%
k = S - delta;
pn = p(neighb_pos);
del_p = [dpdx_p.*fx + (1-fx).*dpdx_p(neighb_pos),dpdy_p.*fx + (1-fx).*dpdy_p(neighb_pos)];%[DPDX(faces1),DPDY(faces1)];
bd = Elements.faces.bound_flag(:,i);
db = faces - centroids;
dn = dot(db,S,2).*S./(S(:,1).^2 + S(:,2).^2);
k(bd==1,:)=0;
neighb_bound = a.neighb_bound(:,i);
pn(bd==1,:)=p_bound(neighb_bound(bd==1));
d(bd==1,:)=dn(bd==1,:);
delta(bd==1,:)=S(bd==1,:);
e = d./vecnorm(d,2,2);



del_p = [del_p(:,1) + ((pn-p)./vecnorm(d,2,2) - dot(del_p,e,2)).*e(:,1),...
          del_p(:,2) + ((pn-p)./vecnorm(d,2,2) - dot(del_p,e,2)).*e(:,2)];

uf = uf + dt/rho * ((pn-p)./vecnorm(d,2,2) - dot(del_p,e,2)).*e(:,1);
vf = vf + dt/rho * ((pn-p)./vecnorm(d,2,2) - dot(del_p,e,2)).*e(:,2);
      
      
uf(bd==1) = u_bound(neighb_bound(bd==1));
vf(bd==1) = v_bound(neighb_bound(bd==1));
      
rhs = rhs + rho/dt *(uf.*nx.*ds + vf.*ny.*ds);

k_delp = k_delp + dot(k,del_p,2);
pn_delta_d = pn_delta_d + pn.*vecnorm(delta,2,2)./vecnorm(d,2,2);
delta_d = delta_d + vecnorm(delta,2,2)./vecnorm(d,2,2);



j=j+2;      
      
end

p_centre = (k_delp - rhs + pn_delta_d)./delta_d;

end
