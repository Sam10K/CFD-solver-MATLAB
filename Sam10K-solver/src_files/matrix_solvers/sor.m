function[x_new,x_init_res,diff,iter] = sor(A,B,x0,tol,max_iter)


D = sparse(diag(diag(A)));
L = sparse(tril(A,-1));
U = sparse(triu(A,1)); 

x_new=x0;w=1;diff=1;iter=0;
while (diff>tol) && (iter<max_iter)  
    
   x_sor = (D+w*L) \ (w*B - (w*U + (w-1)*D)*x_new);
   diff =  norm((D+w*L)*x_new - (w*B - (w*U + (w-1)*D)*x_new))/norm((D+w*L)*x_new);
   x_new=1*x_sor+0*x_new;
   iter=iter+1;
   if iter==2 || iter==1
       x_init_res = diff;  
   end
    
end

% if iter>=max_iter
%     disp('Max iterations reached')
% end
end