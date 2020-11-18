function[s]=tecwriter(var,Nodes_xy,var_names,e2n,filename,time)

tdata = struct;
tdata.Nvars = size(var,2);
tdata.varnames=var_names;
tdata.FEsurfaces.e2n = e2n;
tdata.FEsurfaces.x = Nodes_xy(:,1);
tdata.FEsurfaces.y = Nodes_xy(:,2);
tdata.FEsurfaces.order = 3;
tdata.FEsurfaces.varloc = 1;
tdata.FEsurfaces.v = var';
tdata.FEsurfaces.solutiontime = time;

s=mat2tecplot(tdata,filename);

end

