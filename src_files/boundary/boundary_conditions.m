function[Boundaries] = boundary_conditions(Boundaries,u,v,p,BC)

names = fieldnames(BC);
loc_bound=[];u_bound=[];v_bound=[];p_bound=[];
ubound_type=[];vbound_type=[];pbound_type=[];bound_face_id=[];
for i=1:length(names)
    
    if strcmp(names{i},'internal')
        continue
    end
    
    type = BC.(names{i}).U.type;
    value = BC.(names{i}).U.value;
    if strcmp(type,'fixed')
        U = value*ones(size(Boundaries.(names{i}).faces.cells));
        Boundaries.(names{i}).u = value*ones(size(Boundaries.(names{i}).faces.cells));
        utype = ones(size(Boundaries.(names{i}).faces.cells));
    elseif strcmp(type,'zeroGradient')
        U = u(Boundaries.(names{i}).faces.cells);
        Boundaries.(names{i}).u = u(Boundaries.(names{i}).faces.cells);
        utype = 2*ones(size(Boundaries.(names{i}).faces.cells));
    elseif strcmp(type,'wall')
        U = value*ones(size(Boundaries.(names{i}).faces.cells));
        Boundaries.(names{i}).u = value*ones(size(Boundaries.(names{i}).faces.cells));
        utype = 3*ones(size(Boundaries.(names{i}).faces.cells));
    end
    
    type = BC.(names{i}).V.type;
    value = BC.(names{i}).V.value;
    if strcmp(type,'fixed')
        V = value*ones(size(Boundaries.(names{i}).faces.cells));
        Boundaries.(names{i}).v = value*ones(size(Boundaries.(names{i}).faces.cells));
        vtype = ones(size(Boundaries.(names{i}).faces.cells));
    elseif strcmp(type,'zeroGradient')
        V = v(Boundaries.(names{i}).faces.cells);
        Boundaries.(names{i}).v = v(Boundaries.(names{i}).faces.cells);
        vtype = 2*ones(size(Boundaries.(names{i}).faces.cells));
    elseif strcmp(type,'wall')
        V = value*ones(size(Boundaries.(names{i}).faces.cells));
        Boundaries.(names{i}).v = value*ones(size(Boundaries.(names{i}).faces.cells));
        vtype = 3*ones(size(Boundaries.(names{i}).faces.cells));
    end

    type = BC.(names{i}).P.type;
    value = BC.(names{i}).P.value;
    if strcmp(type,'fixed')
        P = value*ones(size(Boundaries.(names{i}).faces.cells));
        Boundaries.(names{i}).p = value*ones(size(Boundaries.(names{i}).faces.cells));
        ptype = ones(size(Boundaries.(names{i}).faces.cells));
    elseif strcmp(type,'zeroGradient')
        P = p(Boundaries.(names{i}).faces.cells);
        Boundaries.(names{i}).p = p(Boundaries.(names{i}).faces.cells);
        ptype = 2*ones(size(Boundaries.(names{i}).faces.cells));
    elseif strcmp(type,'wall')
        P = p(Boundaries.(names{i}).faces.cells);
        Boundaries.(names{i}).p = p(Boundaries.(names{i}).faces.cells);
        ptype = 3*ones(size(Boundaries.(names{i}).faces.cells));
    end
    
    loc = Boundaries.(names{i}).faces.mid;
    id = Boundaries.(names{i}).faces.id;
    loc_bound = [loc_bound;loc];
    bound_face_id = [bound_face_id;id];
    u_bound = [u_bound;U];
    v_bound = [v_bound;V];
    p_bound = [p_bound;P];
    ubound_type = [ubound_type;utype];
    vbound_type = [vbound_type;vtype];
    pbound_type = [pbound_type;ptype];
    
end

Boundaries.u_bound = u_bound;
Boundaries.v_bound = v_bound;
Boundaries.p_bound = p_bound;
Boundaries.loc_bound = loc_bound;
Boundaries.bound_face_id = bound_face_id;

Boundaries.ubound_type = ubound_type;
Boundaries.vbound_type = vbound_type;
Boundaries.pbound_type = pbound_type;


end