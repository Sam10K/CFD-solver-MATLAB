function[ugrad,vgrad,pgrad] = gradient_loop(Elements,Boundaries,u,v,p,Schemes)

grad_scheme = Schemes.grad_scheme;

if strcmp(grad_scheme,'LSQ')
    [ugrad,vgrad,pgrad] = gradient_loop_lsq(Elements,Boundaries,u,v,p);
elseif strcmp(grad_scheme,'GC')
    [ugrad,vgrad,pgrad] = gradient_loop_gauss_cell(Elements,Boundaries,u,v,p);
elseif strcmp(grad_scheme,'GN')
    [ugrad,vgrad,pgrad] = gradient_loop_gauss_node(Elements,Boundaries,u,v,p);
end

end