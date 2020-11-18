function[x,init_res,res,iter] = matrix_solve(A,B,x0,solver_params)

solver = solver_params.solver;
max_iter = solver_params.max_iter;
tol = solver_params.tol;

if strcmp(solver,'SOR')
    [x,init_res,res,iter] = sor(A,B,x0,tol,max_iter);
elseif strcmp(solver,'JAC')
    [x,init_res,res,iter] = jacobi(A,B,x0,tol,max_iter);
elseif strcmp(solver,'ILU')
    [x,init_res,res,iter] = ilu_solve(A,B,x0,tol,max_iter);
end



end