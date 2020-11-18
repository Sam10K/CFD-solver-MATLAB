function[] = steady_solver(Solver_params,SIMPLE_params,Schemes)

global Elements;global Boundaries;global Properties;
global u;global v;global p; global BC;


ustar = u; vstar = v;

if Solver_params.residuals=='T'
    figure(1)
    u_res = animatedline('Color','r');
    v_res = animatedline('Color','b');
    p_res = animatedline('Color','k');
    title('Residuals')
    legend('U','V','P')
    xlabel('Time')
    ylabel('Residual')
end

[Elements] = face_interp(Elements,Boundaries,u,v);

iter_outer=0; alpha_u = SIMPLE_params.alpha_u; alpha_p = SIMPLE_params.alpha_p;
while iter_outer<Solver_params.T
    
    pdash = zeros(size(p));

   if Solver_params.flowvis=='T'
        figure(2)
        quiver(Elements.centroid(:,1),Elements.centroid(:,2),u,v);
        pause(0.1);
   end
   
   [~,~,grad_p] = gradient_loop(Elements,Boundaries,ustar,vstar,p,Schemes);
   
   for i=1:Schemes.u_matrix_solver.nsweeps

       [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,p,BC);

       [Ap_x,H_x,~,~] = H_Ap_assemble(Boundaries,Elements,Properties,ustar,vstar,p,Schemes);
       
       [A,B] = eq_relax(Ap_x,H_x - Elements.volume.*grad_p(:,1),ustar,SIMPLE_params.alpha_u);
       
       [ustar_new,u_init_res,diff,iter] = matrix_solve(A,B,ustar,Schemes.u_matrix_solver);
       
       ustar = ustar_new;
       
   end
   
   for i=1:Schemes.v_matrix_solver.nsweeps

       [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,p,BC);

       [~,~,Ap_y,H_y] = H_Ap_assemble(Boundaries,Elements,Properties,ustar,vstar,p,Schemes);
       
       [A,B] = eq_relax(Ap_y,H_y - Elements.volume.*grad_p(:,2),vstar,SIMPLE_params.alpha_u);

       [vstar_new,v_init_res,diff,iter] = matrix_solve(A,B,vstar,Schemes.v_matrix_solver);

       vstar = vstar_new;      
       
   end

   
   [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,p,BC);
   [Ap_x,H_x,Ap_y,H_y] = H_Ap_assemble(Boundaries,Elements,Properties,ustar,vstar,p,Schemes);
   Ap = [full(diag(Ap_x))/alpha_u./Elements.volume,full(diag(Ap_y))/alpha_u./Elements.volume];
   [Elements] = rhie_chow(Elements,Boundaries,ustar,vstar,p,Ap,alpha_u,Schemes);
   
   
   for i=1:Schemes.p_matrix_solver.nsweeps  
       
       [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,pdash,BC);

       pdash_new = pressure_poisson_fract_loop(Boundaries,Elements,ustar,vstar,pdash,Ap,Schemes.p_matrix_solver,Schemes);

       diff_pdash = norm(pdash_new-pdash)/norm(pdash_new);

       pdash = pdash_new;
       
   end
   
   [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,pdash,BC);

   [~,~,grad_p] = gradient_loop(Elements,Boundaries,ustar,vstar,pdash,Schemes);
   
   udash = - grad_p(:,1)./Ap(:,1);
   
   vdash = - grad_p(:,2)./Ap(:,2);
   
   ustar = ustar + udash;
   
   vstar = vstar + vdash;
   
   p = p + alpha_p*pdash;
   
   iter_outer = iter_outer + 1;
  
   u = ustar; v = vstar;

   [Boundaries] = boundary_conditions(Boundaries,u,v,p,BC);
   
   [Elements] = rhie_chow(Elements,Boundaries,u,v,p,Ap,alpha_u,Schemes);
   
   [Cont] = continuity(Elements);

   cont = sum(abs(Cont).*Elements.volume)/sum(Elements.volume);
   
   ures = norm(udash)/norm(u);
   vres = norm(vdash)/norm(v);
   pres = norm(pdash)/norm(p);
   
   if isnan(ures)
       break
   end
   
   utol = SIMPLE_params.u_tol;
   vtol = SIMPLE_params.v_tol;
   ptol = SIMPLE_params.p_tol;
      
   if Solver_params.residuals=='T'
       addpoints(u_res,iter_outer,ures);
       addpoints(v_res,iter_outer,vres);
       addpoints(p_res,iter_outer,pres);
       drawnow limitrate
   end
   
   fprintf('U_res = %g, V_res = %g, P_res = %g\n',ures,vres,pres);

   fprintf('T = %g, continuity = %g\n\n',iter_outer,cont);
   
   if Solver_params.write_files=='T' && iter_outer>=Solver_params.write_start && ...
      iter_outer<=Solver_params.write_end && ...
      mod(iter_outer-Solver_params.write_start,Solver_params.write_freq)==0
          
          fprintf('Writing files........\n\n')
          e2n = [Elements.faces.nodes];
          var_names={'x','y','z','u','v','p'};
          var = [u,v,p];
          Nodes_xy = Elements.Nodes.coord;
          mkdir(string(iter_outer));
          output_loc = strcat(string(iter_outer),'/');
          tecname = strcat(output_loc,sprintf('%g.plt',iter_outer));
          s = tecwriter(var,Nodes_xy,var_names,e2n,tecname,iter_outer);
          
   end
  
   
   if ures<=utol && vres<=vtol && pres<=ptol
       disp('Solution converged!')
       break
   end
   
end

end





% if ~isempty(solver_params.output_loc)
%     
%    e2n = [Elements.faces.nodes];
%    var_names={'x','y','z','u','v','p'};
%    var = [u,v,p];
%    Nodes_xy = Nodes;
%    tecname = strcat(output_loc,sprintf('%g.plt',iter_outer));
%    s = tecwriter(var,Nodes_xy,var_names,e2n,tecname,iter_outer);
% 
% end


