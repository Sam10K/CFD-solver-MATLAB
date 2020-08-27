%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the diffusion term of the Navier-stokes over the
% entire domain, i.e. integral(nudel2(u) dv)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[diff] = diffusion(Boundaries,Elements,Properties,u,v,p)

nu = Properties.nu;

a = Elements.faces;
[~,nd] = size(a.nodes);
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;

centroids = [Elements.centroid];
[Ugrad,Vgrad,~] = gradient_loop(Elements,Boundaries,u,v,p);
dudx_p = Ugrad(:,1);
dudy_p = Ugrad(:,2);
dvdx_p = Vgrad(:,1);
dvdy_p = Vgrad(:,2);

j=1;
x_diff = zeros(length(centroids),1);
y_diff = x_diff;

for i=1:nd
    
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

x_diff = x_diff + nu*(vecnorm(delta,2,2).*(un-up)./vecnorm(d,2,2) + dot(k,del_u,2));
y_diff = y_diff + nu*(vecnorm(delta,2,2).*(vn-vp)./vecnorm(d,2,2) + dot(k,del_v,2));
      
j=j+2;
     
end

diff = [x_diff,y_diff];
end
