function[F_x,F_y] = forces(obj)

global Elements;global Boundaries;global Properties;
global u;global v;global p;

mu = Properties.mu;

b_cell_id = Boundaries.(obj).faces.cells;
b_face_id = Boundaries.(obj).faces.id;
face_id = Elements.faces.id;
bd = ismember(face_id,b_face_id);
bd = bd(b_cell_id,:);

nx = Elements.faces.normal(b_cell_id,[1:2:end-1]);
ny = Elements.faces.normal(b_cell_id,[2:2:end]);
area = Elements.faces.area(b_cell_id,:);

faces = Boundaries.(obj).faces.mid;
centroids = Elements.centroid(b_cell_id);

bd_col = reshape(bd',[1,numel(bd)])';
nx_col = reshape(nx',[1,numel(nx)])';
ny_col = reshape(ny',[1,numel(ny)])';
area_col = reshape(area',[1,numel(area)])';
nx = -nx_col(bd_col); % To get outward normal
ny = -ny_col(bd_col); % To get outward normal
ds = area_col(bd_col);

uc = u(b_cell_id);
vc = v(b_cell_id);
pc = p(b_cell_id);

ub = Boundaries.(obj).u;
vb = Boundaries.(obj).v;
pb = Boundaries.(obj).p;

S = [nx.*ds,ny.*ds];
db = (faces - centroids);
dn = dot(db,S,2).*S./(S(:,1).^2 + S(:,2).^2);
e = dn./vecnorm(dn,2,2);


b_grad_u = mu*[(ub-uc)./vecnorm(dn,2,2) .*e(:,1),...
            (ub-uc)./vecnorm(dn,2,2) .*e(:,2)];
        
b_grad_v = mu*[(vb-vc)./vecnorm(dn,2,2) .*e(:,1),...
            (vb-vc)./vecnorm(dn,2,2) .*e(:,2)];

        
tau_x = 2*b_grad_u(:,1).*nx + (b_grad_u(:,2)+b_grad_v(:,1)).*ny;
tau_y = (b_grad_u(:,2)+b_grad_v(:,1)).*nx + 2*b_grad_v(:,2).*ny;

p_x = pb.*nx;
p_y = pb.*ny;


F_x = sum((tau_x-p_x).*ds);
F_y = sum((tau_y-p_y).*ds);


end