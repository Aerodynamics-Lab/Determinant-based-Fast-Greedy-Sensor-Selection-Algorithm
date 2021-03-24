function[isensors]= subsensor_DG_r(U, r, p)

[n,~]=size(U);

THETApp=zeros(p,r);
THETAppTHETApptinv=zeros(p,p);
abs=zeros(1,n);
tmp=zeros(1,r);
sensor=zeros(1,p);

for pp=1:p
    Y=eye(r)-THETApp(1:pp-1,:)'*THETAppTHETApptinv(1:pp-1,1:pp-1)*THETApp(1:pp-1,:);
    for nn=1:n
        tmp(1,:)=U(nn,:);
        abs(nn)=tmp(1,:)*Y*tmp(1,:)';
    end

[~,sensor(pp)]=max(abs);%sensorベクトルのpp番目の要素にabs最大になる要素の番号を入れる

THETApp(pp,:)=U(sensor(pp),:);%THETAppに行を追加

THETAppTHETApptinv(1:pp,1:pp)=inv(THETApp(1:pp,:)*THETApp(1:pp,:)');
end
isensors=sensor';
end