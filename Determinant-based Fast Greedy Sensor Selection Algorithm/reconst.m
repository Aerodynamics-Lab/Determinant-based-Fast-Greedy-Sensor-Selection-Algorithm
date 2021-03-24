function [Xestimate]= reconst(Xorg,H,U, r)
z=H*Xorg;
C=H*U(:,1:r);
Xestimate=pinv(C)*z; %B = pinv(A) returns the Moore-Penrose Pseudoinverse of matrix A.