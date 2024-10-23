%% FUNC TEST
P=@(k1,k2,A,m,n,t) exp(-k1.*t).*(A+k2.*(m.^n).*(k1.^(-(n+1))).*((-1).^n).*(gammainc(-k1.*t,n+1)-gammainc(0,n+1)));
PGamma=@(k1,k2,A,m,n,t) (gammainc(-k1.*t,n+1)-gammainc(0,n+1));
P2=@(k1,k2,A,m,n,t)  (k1.^(-(n+1))).*((-1).^n).*(gammainc(-k1.*t,n+1)-gammainc(0,n+1));
%k-=k1, k+=k2
A=0;
n=7;
m=7;
k1=0.45;
k2F=@(m,n) 0.45./(m.^n);

SETx=0:0.01:1;%time /min
figure
for n=5:1:10
    hold on
    plot(SETx,P2(k1,k2F(m,n),A,m,n,SETx))
end
% legend(string(5:1:10))