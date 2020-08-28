%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate the gradient of the three flow variables over the
% entire domain. Calculated using Green-Gauss cell based method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[ugrad,vgrad,pgrad] = gradient(Elements,Boundaries,u,v,p)

a = Elements.faces;
[~,nd] = size(a.nodes);
centroids = Elements.centroid;
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;
p_bound = Boundaries.p_bound;
vol = Elements.volume;

j=1;
u_x_grad = zeros(length(centroids),1);
v_x_grad = u_x_grad;
p_x_grad = u_x_grad;

u_y_grad = u_x_grad;
v_y_grad = u_x_grad;
p_y_grad = u_x_grad;


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
    uf = 0.5*(u+u(neighb_pos));
    vf = 0.5*(v+v(neighb_pos));
    pf = 0.5*(p+p(neighb_pos));
%     uf = fx.*u + (1-fx).*u(neighb_pos);
%     vf = fx.*v + (1-fx).*v(neighb_pos);
%     pf = fx.*p + (1-fx).*p(neighb_pos);
    bd = a.bound_flag(:,i);
    uf(bd==1) = u_bound(neighb_bound(bd==1));
    vf(bd==1) = v_bound(neighb_bound(bd==1));
    pf(bd==1) = p_bound(neighb_bound(bd==1));

    u_x_grad = u_x_grad + (uf.*ds.*nx)./vol; 
    u_y_grad = u_y_grad + (uf.*ds.*ny)./vol; 

    v_x_grad = v_x_grad + (vf.*ds.*nx)./vol; 
    v_y_grad = v_y_grad + (vf.*ds.*ny)./vol; 

    p_x_grad = p_x_grad + (pf.*ds.*nx)./vol; 
    p_y_grad = p_y_grad + (pf.*ds.*ny)./vol; 

    j=j+2;

end

ugrad = [u_x_grad,u_y_grad];
vgrad = [v_x_grad,v_y_grad];
pgrad = [p_x_grad,p_y_grad];

diff=1;tol=1e-5;iter=0;
while diff>tol && iter<=10
j=1;
u_x_grad1 = zeros(length(centroids),1);
v_x_grad1 = u_x_grad1;
p_x_grad1 = u_x_grad1;

u_y_grad1 = u_x_grad1;
v_y_grad1 = u_x_grad1;
p_y_grad1 = u_x_grad1;


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
    r = faces - 0.5*(centroids+neighb);
    uf = 0.5*(u+u(neighb_pos)) + 0.5*dot(ugrad+ugrad(neighb_pos),r,2);
    vf = 0.5*(v+v(neighb_pos)) + 0.5*dot(vgrad+vgrad(neighb_pos),r,2);
    pf = 0.5*(p+p(neighb_pos)) + 0.5*dot(pgrad+pgrad(neighb_pos),r,2);
%     uf = fx.*u + (1-fx).*u(neighb_pos);
%     vf = fx.*v + (1-fx).*v(neighb_pos);
%     pf = fx.*p + (1-fx).*p(neighb_pos);
    bd = a.bound_flag(:,i);
    uf(bd==1) = u_bound(neighb_bound(bd==1));
    vf(bd==1) = v_bound(neighb_bound(bd==1));
    pf(bd==1) = p_bound(neighb_bound(bd==1));

    u_x_grad1 = u_x_grad1 + (uf.*ds.*nx)./vol; 
    u_y_grad1 = u_y_grad1 + (uf.*ds.*ny)./vol; 

    v_x_grad1 = v_x_grad1 + (vf.*ds.*nx)./vol; 
    v_y_grad1 = v_y_grad1 + (vf.*ds.*ny)./vol; 

    p_x_grad1 = p_x_grad1 + (pf.*ds.*nx)./vol; 
    p_y_grad1 = p_y_grad1 + (pf.*ds.*ny)./vol; 

    j=j+2;

end

ugrad_new = [u_x_grad1,u_y_grad1];
vgrad_new = [v_x_grad1,v_y_grad1];
pgrad_new = [p_x_grad1,p_y_grad1];

diff = norm([norm(ugrad_new(:,1)-ugrad(:,1)),norm(ugrad_new(:,2)-ugrad(:,2)),...
            norm(vgrad_new(:,1)-vgrad(:,1)),norm(vgrad_new(:,2)-vgrad(:,2)),...
            norm(pgrad_new(:,1)-pgrad(:,1)),norm(pgrad_new(:,2)-pgrad(:,2))]);

iter=iter+1;
ugrad = ugrad_new;
vgrad = vgrad_new;
pgrad = pgrad_new;

end

end