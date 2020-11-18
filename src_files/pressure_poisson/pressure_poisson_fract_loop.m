function[p_centre] = pressure_poisson_fract_loop(Boundaries,Elements,u,v,p,Ap,solver_params,Schemes)

a = Elements.faces;
[~,nd] = size(a.nodes);
p_bound = Boundaries.p_bound;

centroids = [Elements.centroid];
[~,~,Pgrad] = gradient_loop(Elements,Boundaries,u,v,p,Schemes);

j=1;

rhs = zeros(length(centroids),1);
A = sparse(size(p,1),size(p,1)); B = rhs;

for i=1:nd
    

nx = a.normal(:,j);
ny = a.normal(:,j+1);
ds = a.area(:,i);
faces = a.mid(:,[j,j+1]);
S = [ds.*nx,ds.*ny];
neighb_pos = a.neighb(:,i);
neighb = centroids(neighb_pos,:);
bd = logical(Elements.faces.bound_flag(:,i));
pbound_type = a.pbound_type(:,i);
bdzg = pbound_type==2;
bdw = pbound_type==3;
neighb_bound = a.neighb_bound(:,i);
fx = vecnorm(faces-neighb,2,2)./vecnorm(centroids-neighb,2,2);
fx(abs(fx)==Inf)=0;
fx(isnan(fx)) = 0;

Ap_f = [fx.*Ap(:,1) + (1-fx).*Ap(neighb_pos,1),...
        fx.*Ap(:,2) + (1-fx).*Ap(neighb_pos,2)];
Ap_f(bd,:) = Ap(bd,:);

S = [1./Ap_f(:,1) .*S(:,1),1./Ap_f(:,2) .*S(:,2)];

d = (neighb-centroids);
delta = d./dot(d,S,2) .* vecnorm(S,2,2).^2;
k = S - delta;
pn = p(neighb_pos);
db = (faces - centroids);
dn = dot(db,S,2).*S./(S(:,1).^2 + S(:,2).^2);
k(bd,:)=0;
pn(bd,:)=p_bound(neighb_bound(bd));
d(bd,:)=dn(bd,:);
delta(bd,:)=S(bd,:);

del_p = grad_interp(faces,centroids,neighb,neighb_pos,p,Pgrad,bd,S,neighb_bound,p_bound);
del_p(bdzg,:) = 0;
del_p(bdw,:) = 0;


ap = -vecnorm(delta,2,2)./vecnorm(d,2,2);
an = vecnorm(delta,2,2)./vecnorm(d,2,2);
b = dot(k,del_p,2);

ap(bd) = -vecnorm(delta(bd,:),2,2)./vecnorm(d(bd,:),2,2);
an(bd)=0;
b(bd)=pn(bd).*vecnorm(delta(bd,:),2,2)./vecnorm(d(bd,:),2,2);

ap(bdzg)=0;
an(bdzg)=0;
b(bdzg)=0;

ap(bdw)=0;
an(bdw)=0;
b(bdw)=0;

mat = sparse([1:size(p,1)]',[1:size(p,1)]',ap,size(p,1),size(p,1)) + ...
      sparse([1:size(p,1)]',neighb_pos,an,size(p,1),size(p,1));

uf = Elements.faces.uf(:,i);
vf = Elements.faces.vf(:,i);
rhs = (uf.*nx.*ds + vf.*ny.*ds) - b;
rhs(bdw)=0;

A = A + mat;
B = B + rhs;    

j=j+2;      
   
end

[p_centre,p_init_res,diff,iter] = matrix_solve(A,B,p,solver_params);


end
