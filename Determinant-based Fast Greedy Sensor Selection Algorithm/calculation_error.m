function  [Xestimate_rand, Xestimate_convex, Xestimate_QR,Xestimate_DG,Xestimate_QD,...
    Perror_rand,Perror_std_rand,Perror_convex,Perror_std_convex,Perror_QR,Perror_std_QR,Perror_DG,Perror_std_DG,Perror_QD,Perror_std_QD]=...
    calculation_error(v, q, w, r, TT, Xorg, Psi, H_rand, H_convex, H_QR, H_DG, H_QD)
% Rand
[Xestimate_rand]= reconst(Xorg,H_rand,Psi,r);
[Perror_rand,Perror_std_rand] = subcalculation_error(TT,Xorg, Psi, Xestimate_rand,r);

% convex
[Xestimate_convex]= reconst(Xorg,H_convex,Psi,r);
[Perror_convex,Perror_std_convex] =  subcalculation_error(TT,Xorg, Psi, Xestimate_convex,r);
% QR
[Xestimate_QR]= reconst(Xorg,H_QR,Psi,r);
[Perror_QR,Perror_std_QR] =  subcalculation_error(TT,Xorg, Psi, Xestimate_QR,r);
% DG
[Xestimate_DG]= reconst(Xorg,H_DG,Psi,r);
[Perror_DG,Perror_std_DG] =  subcalculation_error(TT,Xorg, Psi, Xestimate_DG,r);
% QD
[Xestimate_QD]= reconst(Xorg,H_QD,Psi,r);
[Perror_QD,Perror_std_QD] =  subcalculation_error(TT,Xorg, Psi, Xestimate_QD,r);
end