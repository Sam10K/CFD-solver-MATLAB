function[SIMPLE_params] = read_SIMPLE_params()

file = './solver_setup/SIMPLE_params';
fid = fopen(file,'r');
data = fgetl(fid);

SIMPLE_params = struct;

while ~feof(fid)
    
    alpha_p = sscanf(data,'alpha_p = %s;');
    if ~isempty(alpha_p)
        alpha_p = split(alpha_p,';');
        SIMPLE_params.alpha_p = str2double(alpha_p(1));
    end
    
    alpha_u = sscanf(data,'alpha_u = %s;');
    if ~isempty(alpha_u)
        alpha_u = split(alpha_u,';');
        SIMPLE_params.alpha_u = str2double(alpha_u(1));
    end
    
    u_tol = sscanf(data,'u_tol = %s;');
    if ~isempty(u_tol)
        u_tol = split(u_tol,';');
        SIMPLE_params.u_tol = str2double(u_tol(1));
    end
    
    v_tol = sscanf(data,'v_tol = %s;');
    if ~isempty(v_tol)
        v_tol = split(v_tol,';');
        SIMPLE_params.v_tol = str2double(v_tol(1));
    end
    
    p_tol = sscanf(data,'p_tol = %s;');
    if ~isempty(p_tol)
        p_tol = split(p_tol,';');
        SIMPLE_params.p_tol = str2double(p_tol(1));
    end
    
    Max_iter = sscanf(data,'Max_iter = %s;');
    if ~isempty(Max_iter)
        Max_iter = split(Max_iter,';');
        SIMPLE_params.Max_iter = str2double(Max_iter(1));
    end
    
    
    
    data = fgetl(fid);
    
end

fclose(fid);

   
end