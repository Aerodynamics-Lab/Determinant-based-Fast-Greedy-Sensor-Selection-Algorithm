    function[isensors]= sensor_GR(U, r, p, s, Xorg)
    if r==p*s
    else
        U=U*U';
    end
        [n,~]=size(U);

        C=[];
        for i=1:s
        C(:,:,i)=U((1+(i-1)*n/s):(i*n/s) ,1:end); 
        end
       % C(:,:,1)=U(1:n/s    ,1:r);
        %C(:,:,2)=U((n/s+1):n,1:r);
        for pp=1:p
            for nn=1:n/s
                Cabs(nn)=1;%C(nn,:,1)*(C(nn,:,1))';
                for ss=1:s
                    Ctmp(ss,:)=C(nn,:,ss);
                end
                for ss=1:s
                     Cabs(nn)=Cabs(nn)*Ctmp(ss,:)*Ctmp(ss,:)';
                     CC=Ctmp(ss,:)/norm(Ctmp(ss,:));
                    for st=ss+1:s
                        Ctmp(st,:)=(Ctmp(st,:))-(CC*Ctmp(st,:)')*(CC);
                    end
                end
            end
            [~,p(pp)]=max(Cabs);
    %         isensors=p;
            for ss=1:s
                CC=C(p(pp),:,ss);
                CC=CC/norm(CC);
                for st=1:s
                C(:,:,st)=C(:,:,st)-C(:,:,st)*CC'*CC;
                end
            end
       end
       isensors=p';

    end