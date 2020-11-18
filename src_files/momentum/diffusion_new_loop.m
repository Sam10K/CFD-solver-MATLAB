function[diff] = diffusion_new_loop(Boundaries,Elements,Properties,u,v,p,Schemes)

nu = Properties.nu;

a = Elements.faces;
[~,nd] = size(a.nodes);
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;

centroids = [Elements.centroid];
[ugrad,vgrad,~] = gradient_loop(Elements,Boundaries,u,v,p,Schemes);

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
d = (neighb-centroids);
delta = d./dot(d,S,2) .* vecnorm(S,2,2).^2;%d.*vecnorm(S,2,2)./vecnorm(d,2,2);
bd = Elements.faces.bound_flag(:,i);
neighb_bound = a.neighb_bound(:,i);
k = S - delta;
un = u(neighb_pos);
vn = v(neighb_pos);
up = u;
vp = v;

del_u = grad_interp(faces,centroids,neighb,neighb_pos,u,ugrad,bd,S,neighb_bound,u_bound);
del_v = grad_interp(faces,centroids,neighb,neighb_pos,v,vgrad,bd,S,neighb_bound,v_bound);

Del_u_velocity = [del_u(:,1),del_v(:,1)];% First row of Del_u tensor
Del_v_velocity = [del_u(:,2),del_v(:,2)];% Second row of Del_u tensor

db = (faces - centroids);
dn = dot(db,S,2).*S./(S(:,1).^2 + S(:,2).^2);
k(bd==1,:)=0;
un(bd==1,:)=  u_bound(neighb_bound(bd==1));
vn(bd==1,:)=  v_bound(neighb_bound(bd==1));
d(bd==1,:)=dn(bd==1,:);
delta(bd==1,:)=S(bd==1,:);

x_diff = x_diff + nu*(vecnorm(delta,2,2).*(un-up)./vecnorm(d,2,2) + dot(Del_u_velocity,k,2) + dot(del_u,S,2));
y_diff = y_diff + nu*(vecnorm(delta,2,2).*(vn-vp)./vecnorm(d,2,2) + dot(Del_v_velocity,k,2) + dot(del_v,S,2));
      
j=j+2;
     
end

diff = [x_diff,y_diff];
end
