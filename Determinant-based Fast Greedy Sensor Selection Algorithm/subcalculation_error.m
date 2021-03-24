function[Perror_ave,Perror_std] = subcalculation_error(TT,Xorg, Psi, Xestimate,r)

for ii=1:TT
    x=Xorg(:,ii);
    xes=Psi(:,1:r)*Xestimate(:,ii);
    Perror(ii) = norm(xes-x)/norm(x);
end
Perror_ave=mean(Perror);
Perror_std=std(Perror);
end
