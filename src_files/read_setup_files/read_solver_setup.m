function[solver_setup] = read_solver_setup()

file = './solver_setup/solver_setup';
fid = fopen(file,'r');
data = fgetl(fid);

solver_setup = struct;

while ~feof(fid)
    
    steady_flag = sscanf(data,'steady = %c;');
    if ~isempty(steady_flag)
        solver_setup.steady_flag = steady_flag;
    end
    
    mu = sscanf(data,'mu = %s;');
    if ~isempty(mu)
        mu = split(mu,';');
        solver_setup.mu = str2double(mu(1));
    end
    
    rho = sscanf(data,'rho = %s;');
    if ~isempty(rho)
        rho = split(rho,';');
        solver_setup.rho = str2double(rho(1));
    end
    
    T = sscanf(data,'T = %g;');
    if ~isempty(T)
        solver_setup.T = T;
    end
    
    CFL = sscanf(data,'CFL = %g;');
    if ~isempty(CFL)
        solver_setup.CFL = CFL;
    end
    
    flowvis = sscanf(data,'flowvis = %s;');
    if ~isempty(flowvis)
        flowvis = split(flowvis,';');
        solver_setup.flowvis = flowvis{1};
    end
    
    residuals = sscanf(data,'residuals = %s;');
    if ~isempty(residuals)
        residuals = split(residuals,';');
        solver_setup.residuals = residuals{1};
    end
    
    write_files = sscanf(data,'write_files = %c;');
    if ~isempty(write_files)
        solver_setup.write_files = write_files;
    end
    
    write_start = sscanf(data,'write_start = %g;');
    if ~isempty(write_start)
        solver_setup.write_start = write_start;
    end
    
    write_end = sscanf(data,'write_end = %g;');
    if ~isempty(write_end)
        solver_setup.write_end = write_end;
    end
    
    write_freq = sscanf(data,'write_freq = %g;');
    if ~isempty(write_freq)
        solver_setup.write_freq = write_freq;
    end
    
    data = fgetl(fid);
    
end

fclose(fid);

solver_setup.nu = solver_setup.mu/solver_setup.rho;

   
end