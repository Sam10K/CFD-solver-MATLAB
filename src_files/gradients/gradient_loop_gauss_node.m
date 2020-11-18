function[ugrad,vgrad,pgrad] = gradient_loop_gauss_node(Elements,Boundaries,u,v,p)

a = Elements.faces;
[~,nd] = size(a.nodes);
centroids = Elements.centroid;
u_bound = Boundaries.u_bound;
v_bound = Boundaries.v_bound;
p_bound = Boundaries.p_bound;
vol = Elements.volume;
Node_coord = Elements.Nodes.coord;
Node_cells = Elements.Nodes.cells;


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
    neighb_bound = a.neighb_bound(:,i);
    
    face_nodes = Elements.faces.face_nodes(:,[j,j+1]);
    n1 = face_nodes(:,1);
    n2 = face_nodes(:,2);
    r_n1 = Node_coord(n1,:);
    r_n2 = Node_coord(n2,:);
    cells_n1 = Node_cells(n1); % In cell {} format
    cells_n2 = Node_cells(n2); % In cell {} format
    
    u_n1 = zeros(size(u));v_n1=u_n1;p_n1=u_n1;
    u_n2 = zeros(size(u));v_n2=u_n1;p_n2=u_n1;
    
    for k=1:length(u)
        u_cell_n1 = u(cells_n1{k});
        u_cell_n2 = u(cells_n2{k});
        
        v_cell_n1 = v(cells_n1{k});
        v_cell_n2 = v(cells_n2{k});
        
        p_cell_n1 = p(cells_n1{k});
        p_cell_n2 = p(cells_n2{k});
        
        r_cell_n1 = centroids(cells_n1{k},:);
        r_cell_n2 = centroids(cells_n2{k},:);
        
        dist_n1 = vecnorm(r_n1(k,:)-r_cell_n1,2,2);
        dist_n2 = vecnorm(r_n2(k,:)-r_cell_n2,2,2);
        
        u_n1(k,1) = sum(u_cell_n1./dist_n1)/sum(1./dist_n1);
        u_n2(k,1) = sum(u_cell_n2./dist_n2)/sum(1./dist_n2);
        
        v_n1(k,1) = sum(v_cell_n1./dist_n1)/sum(1./dist_n1);
        v_n2(k,1) = sum(v_cell_n2./dist_n2)/sum(1./dist_n2);
        
        p_n1(k,1) = sum(p_cell_n1./dist_n1)/sum(1./dist_n1);
        p_n2(k,1) = sum(p_cell_n2./dist_n2)/sum(1./dist_n2);
        
        
    end  
    
  
    uf = 0.5*(u_n1+u_n2);
    vf = 0.5*(v_n1+v_n2);
    pf = 0.5*(p_n1+p_n2);

    bd = logical(a.bound_flag(:,i));
    uf(bd) = u_bound(neighb_bound(bd));
    vf(bd) = v_bound(neighb_bound(bd));
    pf(bd) = p_bound(neighb_bound(bd));

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

end