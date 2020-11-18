clc
clear
close

addpath(genpath('./src_files/'))

% function[]=NS_solver()

global Elements;global Boundaries;global Properties; global BC;
global u;global v;global p;

Solver_params = read_solver_setup();
SIMPLE_params = read_SIMPLE_params();
Schemes = read_Schemes();

mesh_file = dir('./mesh_file/');
mesh_file = (strcat(mesh_file(3).folder,'/',mesh_file(3).name));
[Elements,Boundaries,Nodes] = mshread_fluent(mesh_file);

Properties = struct('mu',Solver_params.mu,'rho',Solver_params.rho,'nu',Solver_params.nu);

BC = bc_read();


[u,v,p,Boundaries,Elements] = initialize(BC,Elements,Boundaries);

if Solver_params.steady_flag == 'T'
    
    steady_solver(Solver_params,SIMPLE_params,Schemes);
    
else
    
    unsteady_solver(Solver_params,SIMPLE_params,Schemes);
    
end
% end
    
