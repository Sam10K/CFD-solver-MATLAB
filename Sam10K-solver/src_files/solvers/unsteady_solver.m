function[] = unsteady_solver(Solver_params,SIMPLE_params,Schemes)

global Elements;global Boundaries;global Properties;
global u;global v;global p; global BC;


CFL_max = Solver_params.CFL;
dt = min(CFL_max *mean(Elements.volume./Elements.faces.area,2)/max(u));

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

t=0;
iter_dummy = 0; % Dummy iter counter for writing
while t<=Solver_params.T
    
    conv_n = convection_loop(Boundaries,Elements,u,v,p,Schemes);
    diff_n = diffusion_new_loop(Boundaries,Elements,Properties,u,v,p,Schemes);
    
    fprintf('Time = %g\n\n',t);

    iter_outer=0; alpha_u = SIMPLE_params.alpha_u; alpha_p = SIMPLE_params.alpha_p;
    while iter_outer<SIMPLE_params.Max_iter

        pdash = zeros(size(p));

       if Solver_params.flowvis=='T'
            figure(2)
            quiver(Elements.centroid(:,1),Elements.centroid(:,2),u,v);
            pause(0.1);
       end

       [~,~,grad_p] = gradient_loop(Elements,Boundaries,ustar,vstar,p,Schemes);

       for i=1:Schemes.u_matrix_solver.nsweeps

           [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,p,BC);

           [Ap_x,H_x,~,~] = H_Ap_assemble(Boundaries,Elements,Properties,ustar,vstar,p,Schemes,dt,u,v,conv_n,diff_n);

           [A,B] = eq_relax(Ap_x,H_x - Elements.volume.*grad_p(:,1),ustar,SIMPLE_params.alpha_u);

           [ustar_new,u_init_res,diff,iter] = matrix_solve(A,B,ustar,Schemes.u_matrix_solver);

           ustar = ustar_new;

       end

       for i=1:Schemes.v_matrix_solver.nsweeps

           [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,p,BC);

           [~,~,Ap_y,H_y] = H_Ap_assemble(Boundaries,Elements,Properties,ustar,vstar,p,Schemes,dt,u,v,conv_n,diff_n);

           [A,B] = eq_relax(Ap_y,H_y - Elements.volume.*grad_p(:,2),vstar,SIMPLE_params.alpha_u);

           [vstar_new,v_init_res,diff,iter] = matrix_solve(A,B,vstar,Schemes.v_matrix_solver);

           vstar = vstar_new;      

       end


       [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,p,BC);
       [Ap_x,H_x,Ap_y,H_y] = H_Ap_assemble(Boundaries,Elements,Properties,ustar,vstar,p,Schemes,dt,u,v,conv_n,diff_n);
       Ap = [full(diag(Ap_x))/alpha_u./Elements.volume,full(diag(Ap_y))/alpha_u./Elements.volume];
       [Elements] = rhie_chow(Elements,Boundaries,ustar,vstar,p,Ap,alpha_u,Schemes);


       for i=1:Schemes.p_matrix_solver.nsweeps  

           [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,pdash,BC);

           [pdash_new,p_init_res,diff,iter] = pressure_poisson_fract_loop(Boundaries,Elements,ustar,vstar,pdash,Ap,Schemes.p_matrix_solver,Schemes);

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

       [Boundaries] = boundary_conditions(Boundaries,ustar,vstar,p,BC);

       [Elements] = rhie_chow(Elements,Boundaries,ustar,vstar,p,Ap,alpha_u,Schemes);       

       ures = norm(udash)/norm(ustar);
       vres = norm(vdash)/norm(vstar);
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

       fprintf('\tU_res = %g, V_res = %g, P_res = %g\n',ures,vres,pres);

       if ures<=utol && vres<=vtol && pres<=ptol
           disp('Solution converged!')
           break
       end

    end
    
     u = ustar; v = vstar;
     
     [Boundaries] = boundary_conditions(Boundaries,u,v,p,BC);

     [Elements] = rhie_chow(Elements,Boundaries,u,v,p,Ap,alpha_u,Schemes);
     
     [Cont,Dt] = continuity(Elements,CFL_max);

     cont = sum(abs(Cont).*Elements.volume)/sum(Elements.volume);
     
     dt_new = min(Dt);
     
     t=t+dt;
   
     dt = 0.1*dt_new+0.9*dt;
     
     fprintf('\nT = %g, dt_new = %g, continuity = %g\n\n',t,mean(dt),cont);
     
     if Solver_params.write_files=='T' && t>=Solver_params.write_start && ...
        t<=Solver_params.write_end
    
          
          if mod(iter_dummy,Solver_params.write_freq)==0
              
              fprintf('Writing files........\n\n')
              e2n = [Elements.faces.nodes];
              var_names={'x','y','z','u','v','p'};
              var = [u,v,p];
              Nodes_xy = Elements.Nodes.coord;
              mkdir(string(t));
              output_loc = strcat(string(t),'/');
              tecname = strcat(output_loc,sprintf('%g.plt',t));
              s = tecwriter(var,Nodes_xy,var_names,e2n,tecname,t);
              
              iter_dummy=0;
              
          end
          
          iter_dummy=iter_dummy+1;
          
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


