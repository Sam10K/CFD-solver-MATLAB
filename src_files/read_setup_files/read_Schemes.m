function[Schemes] = read_Schemes()

file = './solver_setup/Schemes';
fid = fopen(file,'r');
data = fgetl(fid);

Schemes = struct;

while ~feof(fid)
    
    div_scheme = sscanf(data,'div_scheme = %s;');
    if ~isempty(div_scheme)
        div_scheme = split(div_scheme,';');
        Schemes.div_scheme = div_scheme{1};
    end
    
    grad_scheme = sscanf(data,'grad_scheme = %s;');
    if ~isempty(grad_scheme)
        grad_scheme = split(grad_scheme,';');
        Schemes.grad_scheme = grad_scheme{1};
    end
    
    u_matrix_solver = sscanf(data,'u_matrix_solver = %s');
    if ~isempty(u_matrix_solver)
        fgetl(fid);curr=fgetl(fid);
        max_iter = sscanf(curr,'max_iter = %d;');
        curr=fgetl(fid);
        tol = sscanf(curr,'tol = %g;');
        curr=fgetl(fid);
        nsweeps = sscanf(curr,'nsweeps = %d;');
        Schemes.u_matrix_solver.solver = u_matrix_solver;
        Schemes.u_matrix_solver.max_iter = max_iter;
        Schemes.u_matrix_solver.tol = tol;
        Schemes.u_matrix_solver.nsweeps = nsweeps;       
    end
    
    v_matrix_solver = sscanf(data,'v_matrix_solver = %s');
    if ~isempty(v_matrix_solver)
        fgetl(fid);curr=fgetl(fid);
        max_iter = sscanf(curr,'max_iter = %d;');
        curr=fgetl(fid);
        tol = sscanf(curr,'tol = %g;');
        curr=fgetl(fid);
        nsweeps = sscanf(curr,'nsweeps = %d;');
        Schemes.v_matrix_solver.solver = v_matrix_solver;
        Schemes.v_matrix_solver.max_iter = max_iter;
        Schemes.v_matrix_solver.tol = tol;
        Schemes.v_matrix_solver.nsweeps = nsweeps;       
    end
    
    p_matrix_solver = sscanf(data,'p_matrix_solver = %s');
    if ~isempty(p_matrix_solver)
        fgetl(fid);curr=fgetl(fid);
        max_iter = sscanf(curr,'max_iter = %d;');
        curr=fgetl(fid);
        tol = sscanf(curr,'tol = %g;');
        curr=fgetl(fid);
        nsweeps = sscanf(curr,'nsweeps = %d;');
        Schemes.p_matrix_solver.solver = p_matrix_solver;
        Schemes.p_matrix_solver.max_iter = max_iter;
        Schemes.p_matrix_solver.tol = tol;
        Schemes.p_matrix_solver.nsweeps = nsweeps;       
    end
    
    
    
    data = fgetl(fid);
    
end

fclose(fid);

   
end