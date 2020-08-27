%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to impose boundary conditions to all the three flow variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Boundaries] = boundary_conditions(Boundaries,u,v,p,BC)

names = fieldnames(BC);
loc_bound=[];u_bound=[];v_bound=[];p_bound=[];

for i=1:length(names)
    
    type = BC.(names{i}).U.type;
    value = BC.(names{i}).U.value;
    if strcmp(type,'fixed')
        U = value*ones(size(Boundaries.(names{i}).faces.cells));
    elseif strcmp(type,'zeroGradient')
        U = u(Boundaries.(names{i}).faces.cells);
    end
    
    type = BC.(names{i}).V.type;
    value = BC.(names{i}).V.value;
    if strcmp(type,'fixed')
        V = value*ones(size(Boundaries.(names{i}).faces.cells));
    elseif strcmp(type,'zeroGradient')
        V = v(Boundaries.(names{i}).faces.cells);
    end

    type = BC.(names{i}).P.type;
    value = BC.(names{i}).P.value;
    if strcmp(type,'fixed')
        P = value*ones(size(Boundaries.(names{i}).faces.cells));
    elseif strcmp(type,'zeroGradient')
        P = p(Boundaries.(names{i}).faces.cells);
    end
    
    loc = Boundaries.(names{i}).faces.mid;
    loc_bound = [loc_bound;loc];
    u_bound = [u_bound;U];
    v_bound = [v_bound;V];
    p_bound = [p_bound;P];
    
end

Boundaries.u_bound = u_bound;
Boundaries.v_bound = v_bound;
Boundaries.p_bound = p_bound;
Boundaries.loc_bound = loc_bound;


end