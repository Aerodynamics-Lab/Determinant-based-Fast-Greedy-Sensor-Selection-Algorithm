clear;
n=1000;
r=10;
ave_number=11;
resi=zeros(100,ave_number);
for rr=1:100
    rr
    for ave=2:ave_number
        U= randn(n,r);
        UUT=U*U';
        [Q,R,pivot] = qr(UUT,'vector');
        U_QR=Q(:,1:rr)*R(1:rr,:);
        
        resi(rr,ave)=norm(Q*R-U_QR);
    end
    resi(rr,1)=mean(resi(rr,2:ave_number));
end