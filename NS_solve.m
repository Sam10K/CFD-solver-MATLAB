%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A 2D Navier-Stokes solver for unsteady,laminar,
% incompressible flows using finite-volume method and collocated grid
% arrangement.
% Pressure-velocity coupling implemented using SIMPLE algorithm.
% Spatial discretization - Central differencing
% Temporal discretization - Implicit Crank-Nicholson
% Accepts all-tri mesh as well as all-quad mesh in fluent mesh file format.
% Parameters to set:
% 1) mesh_file : Mesh file location
% 2) Re : Reynolds number based on the parameter 'L'. Used to compute U_ref
% 3) mu : Fluid dynamic viscosity 
% 4) rho : Fluid density
% 5) L  : Geometric length scale. Usually diameter of the cylinder, airfoil
%                                                        chord length,etc
% 6) CFL_max : Maximum CFL number. Time step size varies based on this
% 7) alpha_p : Under relaxation factor for pressure (usually 0.1-0.3)
% 8) alpha_u : Under relaxation factor for momentum eq (usually 0.3-0.7)
% 9) T : Flow end time
% 10) flowvis : Flag set to visualize flow field every time step as a quiver
%                       plot. Set it to 'T' to activate, 'F' to deactivate
% 11) residuals : flag set to visualize residuals every time step
% 
% Set the boundary conditions using the files 'U.bc', 'V.bc', 'P.bc' 
% inside the folder named BC. Check the example boundary condition files.
% Currently supports fixed value and zero gradient boundary condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
addpath('./solver_functions/')
mesh_file = 'meshes/flat_plate.msh'; % Mesh file location
Re = 10000; % Based on size of square
mu = 1e-4;
rho = 1;
nu = mu/rho;
L = 1;
U_ref = Re*nu/L;
CFL_max = 1;
alpha_p = 0.2;
alpha_u = 0.5;
T = 27400; % End time
flowvis = 'T'; % flag set to visualize flow field every time step
residuals = 'F'; % flag set to visualize residuals

disp('Reading mesh. Please wait ..........')

[Elements,Boundaries,Nodes] = mshread(mesh_file);

fprintf('Mesh reading complete. Mesh has a total of %d elements.\nShowing the mesh in triangulated form......\nPress any key to continue......\n',length(Nodes));

patch('Vertices',[Nodes.x,Nodes.y],'Faces',[Elements.faces.nodes],'FaceColor','None')

pause()
close

Properties = struct('mu',mu,'rho',rho,'nu',nu);


dt_cfl = min(CFL_max*mean(Elements.volume./Elements.faces.area,2)/U_ref);
dt = dt_cfl;

u = U_ref*ones(size(Elements.volume));
v = zeros(size(Elements.volume));
p = zeros(size(Elements.volume));

BC = bc_read();
[Boundaries] = boundary_conditions(Boundaries,u,v,p,BC);
centroids = [Elements.centroid];
loc_bound = [Boundaries.loc_bound];
a = Elements.faces;
[~,nd] = size(a.nodes);
Elements.faces.neighb_bound = zeros(length(centroids),nd);
for i=1:length(centroids)
    j=1;
    for k=1:nd
    
       face = a.mid(i,[j,j+1]);
       dist = (face(1)-loc_bound(:,1)).^2 + (face(2)-loc_bound(:,2)).^2;
       [d,pos] = min(dist);
       if d==0
           Elements.faces.neighb_bound(i,k) = pos;
       else
           Elements.faces.neighb_bound(i,k) = 0;
       end
       j=j+2;
    end

end

disp('Boundary conditions succesfully read. Solver started ........')
ustar = u; vstar = v;


if residuals=='T'
    figure(1)
    u_res = animatedline('Color','r');
    v_res = animatedline('Color','b');
    p_res = animatedline('Color','k');
    title('Residuals')
    legend('U','V','P')
    xlabel('Time')
    ylabel('Residual')
end

t=dt;
% Solver routine starts
while t<=T
    
   if flowvis=='T'
        figure(2)
        quiver(Elements.centroid(:,1),Elements.centroid(:,2),u,v);
        pause(0.1);
        
   end
   
   [~,~,grad_p] = gradient(Elements,Boundaries,u,v,p);
   [Boundaries] = boundary_conditions(Boundaries,u,v,p,BC);
   conv_n = convection(Boundaries,Elements,u,v,p);
   diff_n = diffusion(Boundaries,Elements,Properties,u,v,p);    
       
   
   % Pressure Outer iteration
   diff_outer=1;tol_outer=1e-5;iter_outer=0;
   while (diff_outer>tol_outer) && iter_outer<=50 
       
       
       % Momentum inner iteration
       tol=1e-6;diff=1;iter=0;
       while (diff>tol) && iter<500

           [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,p,BC);

           [H,Ap] = H_Ap(Boundaries,Elements,Properties,ustar,vstar,p,conv_n,diff_n,dt,u,v);

           ustar_new = (H(:,1) - grad_p(:,1))./Ap(:,1);

           diff = rms(ustar_new-ustar)/U_ref;


           iter=iter+1;

           if (iter==1) && (iter_outer==0)
               u_init_res = diff;  
           end

           ustar = alpha_u*ustar_new + (1-alpha_u)*ustar;

       end

%        fprintf('U: Initial residual = %g, Final residual = %g, Iter = %d\n',u_init_res,diff,iter);
       if isnan(diff) break; end

       tol=1e-6;diff=1;iter=0;
       while (diff>tol) && iter<500

           [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,p,BC);

           [H,Ap] = H_Ap(Boundaries,Elements,Properties,ustar,vstar,p,conv_n,diff_n,dt,u,v);

           vstar_new = (H(:,2) - grad_p(:,2))./Ap(:,2);

           diff = rms(vstar_new-vstar)/U_ref;


           iter=iter+1;

           if (iter==1) && (iter_outer==0)
               v_init_res = diff;  
           end

           vstar = alpha_u*vstar_new + (1-alpha_u)*vstar;

       end
%        fprintf('V: Initial residual = %g, Final residual = %g, Iter = %d\n',v_init_res,diff,iter);
       if isnan(diff) break; end
       
       [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,p,BC);
  
       p_new = pressure_poisson(Boundaries,Elements,Properties,ustar,vstar,p,dt);

       diff_outer = rms(p-p_new)/(rho*U_ref^2);

       p = alpha_p*p_new + (1-alpha_p)*p;
  
       [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,p,BC);

       [H,Ap] = H_Ap(Boundaries,Elements,Properties,ustar,vstar,p,conv_n,diff_n,dt,u,v);
              
       [~,~,grad_p] = gradient(Elements,Boundaries,ustar,vstar,p);
              
       ustar = (H(:,1) - grad_p(:,1))./Ap(:,1);
       
       vstar = (H(:,2) - grad_p(:,2))./Ap(:,2);
       
       iter_outer = iter_outer + 1;
       
       if iter_outer==1
            p_init_res = diff_outer;  
       end
       
%        fprintf('P: Initial residual = %g, Final residual = %g, Iter = %d\n',p_init_res,diff_outer,iter_outer);
       if isnan(diff_outer) break; end
   
   
   end
   fprintf('U: Initial residual = %g, Final residual = %g\n',u_init_res,diff);
   fprintf('V: Initial residual = %g, Final residual = %g\n',v_init_res,diff);
   fprintf('P: Initial residual = %g, Final residual = %g, Iter = %d\n',p_init_res,diff_outer,iter_outer);
   
   u = ustar; v = vstar;
   

   [Boundaries] = boundary_conditions(Boundaries,u,v,p,BC);

   [Cont,Dt] = continuity(Boundaries,Elements,u,v,CFL_max);  
   
   cont = norm(Cont);
   
   dt_new = min(Dt);
   
   dt = 0.1*dt_new+0.9*dt;
   
%    dt = min(CFL*mean(Elements.volume./Elements.faces.area,2)./vecnorm([u,v],2,2));
   
   fprintf('T = %g, dt = %g, continuity = %g\n\n',t,dt,cont);
   
   t=t+dt;
   
   if residuals=='T'
       addpoints(u_res,t,u_init_res);
       addpoints(v_res,t,v_init_res);
       addpoints(p_res,t,p_init_res);
       drawnow limitrate
   end
   
    
end

