function[A,B] = eq_relax(A,B,x0,alpha)

diag_col = spdiags(A,0);
A_diag = spdiags(diag_col,0,size(x0,1),size(x0,1));
A = A_diag/alpha + A -A_diag;
B = B +(1-alpha)/alpha * A_diag*x0;


end