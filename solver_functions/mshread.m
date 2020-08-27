%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to read the 2D ascii fluent mesh file and create the structures
% Elements,Boundaries,Nodes to be used in the solver 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Elements,Boundaries,Nodes] = mshread(filename)

fid = fopen(filename,'r');

curr = fgetl(fid);
flag=0;i=1;
while ~strcmp(curr,'))')
    
    if strcmp(curr,"(0 ""Node Section"")")
        flag=1;
        curr=fgetl(fid);
        curr=sscanf(fgetl(fid),'(%x (%x %x %x %x %x)');
        n = curr(4)-curr(3) + 1;
        Nodes = zeros(n,2);
        curr=fgetl(fid);
        curr=fgetl(fid);
    end
    
    if flag        
        
        data = sscanf(curr,'%f %f');
        Nodes(i,:)=data;
        i=i+1;
        
    end
    
    curr = fgetl(fid); 
        
    
end

Boundaries = struct();flag=0;i=1;
Boundaries.loc_bound=[];
face_nodes=[];face_cells=[];
while ~feof(fid)
    
    dummy = sscanf(curr,'(12 (%x %x %x %x %x))');
    if ~isempty(dummy) Ncells = dummy(3)-dummy(2)+1;end
        
   if contains(curr,"zone") 
       
       data = split(curr);
       data = split(data(end),'")');
       Boundaries.(data{1}).faces.mid = [];
       Boundaries.(data{1}).faces.cells = [];
       names{i} = data{1};
       Boundaries.names = names;
       i=i+1;
       fgetl(fid);
       curr = fgetl(fid);
       flag=1;
       
   end
   
   if flag
       
       nodes = sscanf(curr,'%x %x %x %x');
       if ~isempty(nodes)
           
           loc = nodes(1:2);
           x = Nodes(loc,1); y = Nodes(loc,2);
           x_mid = 0.5*(x(1)+x(2));
           y_mid = 0.5*(y(1)+y(2));
           Boundaries.(data{1}).faces.mid = [Boundaries.(data{1}).faces.mid;[x_mid,y_mid]];
           Boundaries.(data{1}).faces.cells = [Boundaries.(data{1}).faces.cells;nodes(3)];
%            Boundaries.loc_bound = [Boundaries.loc_bound;[x_mid,y_mid]];
           
           face_nodes = [face_nodes;nodes(1:2)'];
           face_cells = [face_cells;nodes(3:4)'];
       end
       
       
       
   end
   
   curr = fgetl(fid);
   
end

Elements = struct('id',[],'type',[],'patch',[],'centroid',[],'faces',[],'volume',[],'ghosts',[]);
Elements.faces.nodes=[];
Elements.faces.area=[];
Elements.faces.normal=[];
Elements.faces.mid=[];
cells1 = face_cells(:,1);cells2 = face_cells(:,2);
for i=1:Ncells
   
   face1 = find(i==cells1);face2 = find(i==cells2);
   faces = [face1;face2]';
   node1 = face_nodes(face1,:);node2 = face_nodes(face2,:);
   nodes = reshape([node1;node2]',[1,numel([node1;node2]')]);
   dummy = unique(nodes);
   [loc,vol] = convhull(Nodes(dummy,:));
   
   x = Nodes(dummy(loc(1:end-1)),1);
   y = Nodes(dummy(loc(1:end-1)),2);
   out = polygeom(x,y);
   x_centre = out(2);
   y_centre = out(3);

   k=1;
   for j=1:length(faces)
       
      x = Nodes(face_nodes(faces(j),:),1)';
      y = Nodes(face_nodes(faces(j),:),2)';  
      
      x_mid = 0.5*(x(1)+x(2));
      y_mid = 0.5*(y(1)+y(2));
      
      area = sqrt((x(2)-x(1)).^2 + (y(2)-y(1)).^2);
      
      sign_nx = sign(x_mid - x_centre);
      sign_ny = sign(y_mid - y_centre);
      
      m = (y(2)-y(1))/(x(2)-x(1));
      nx = 1/sqrt(1+1/m^2) * sign_nx;
      ny = 1/sqrt(1+m^2) * sign_ny;
      
      cells = face_cells(faces(j),:);
      neighb = cells(cells~=i);
      
      Elements.faces.area(i,j) = area;
      Elements.faces.mid(i,[k,k+1]) = [x_mid,y_mid];
      Elements.faces.normal(i,[k,k+1]) = [nx,ny];
      Elements.faces.neighb(i,j) = neighb + 1*(neighb==0);
      Elements.faces.bound_flag(i,j) = 1*(neighb==0);
      Elements.faces.neighb_bound(i,j) = faces(j);
      k=k+2;
       
   end
   
   Elements.volume(i,1) = vol;
   Elements.faces.nodes(i,:) = dummy(loc(1:end-1));
   Elements.centroid(i,:) = [x_centre,y_centre];
   Elements.faces.id(i,:) = faces;
   
   
end

fclose(fid);


end

