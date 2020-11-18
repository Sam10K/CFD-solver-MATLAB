function[x_new,x_init_res,diff,iter] = ilu_solve(A,B,x0,tol,max_iter)


[L,U] = ilu(A);
R = A-L*U;

x_new=x0;diff=1;iter=0;
while (diff>tol) && (iter<max_iter)  
    
   x_ilu = (A-R)\((A-R)*x_new + (B-A*x_new));
   diff =  norm(x_ilu-x_new)/norm(x_ilu);
   x_new=1*x_ilu+0*x_new;
   iter=iter+1;
   if iter==1
       x_init_res = diff;  
   end
    
end

% if iter>=max_iter
%     disp('Max iterations reached')
% end
end