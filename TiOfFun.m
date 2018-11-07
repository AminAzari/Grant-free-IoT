%IHN
function D=TiOfFun(Np,t,Tb,r1,r2,Nt,Leng_pack,Tperiod,InSyn)
D=zeros(Nt,Np);
if(InSyn==1)
    Dvec1=((Tperiod-Leng_pack)*rand(1,Np));
    Dvec=floor(Dvec1*length(t)/Tb)/( length(t)/Tb);
    D(1,:)=Dvec; 
    DD=(r1+floor(r2*rand(1,Np)));
    for i=2:Nt
        D(i,:)=D(i-1,:)+[DD .*floor(length(t)/Tb*Leng_pack)/(length(t)/Tb)];
    end
    
elseif(InSyn==2)
    Dvec2=(floor(1+floor(Tperiod/Leng_pack)*rand(1,Np)));
    Dvec1=(Dvec2-1)*Leng_pack;
    Dvec=floor(Dvec1*length(t)/Tb)/( length(t)/Tb);
    D(1,:)=Dvec;
    DD= (r1+floor(r2*rand(1,Np)));
    for i=2:Nt
        D(i,:)=D(i-1,:)+[ DD.*floor(length(t)/Tb*Leng_pack)/(length(t)/Tb)];
    end
end
% for i=1:Np
%     for j=2:Nt
%         if(D(j,i)>Tperiod)
%             D(j,i)=D(1,i);
%         end
%     end
% end
