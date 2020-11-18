function[x_new,x_init_res,diff,iter] = jacobi(A,B,x0,tol,max_iter)


D = sparse(diag(diag(A)));
L = sparse(tril(A,-1));
U = sparse(triu(A,1)); 

x_new=x0;diff=1;iter=0;
while (diff>tol) && (iter<max_iter)  
    
   x_jacobi = D \ (B - (L+U)*x_new);
   diff =  norm(x_jacobi-x_new)/norm(x_jacobi);
   x_new=1*x_jacobi+0*x_new;
   iter=iter+1;
   if iter==1
       x_init_res = diff;  
   end
    
end

% if iter>=max_iter
%     disp('Max iterations reached')
% end
end