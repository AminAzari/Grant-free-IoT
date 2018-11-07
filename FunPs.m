function pp=FunPs(Nres,lambda,Tpr)
KK=10000;
NU=poissrnd(lambda*Tpr,1,KK);
for i=1:KK
    if(NU(i)==0 || NU(i)==1)
        p(i)=1;
    else
        p(i)=((Nres-1)/Nres)^(NU(i)-1);
    end
end
pp=mean(p);