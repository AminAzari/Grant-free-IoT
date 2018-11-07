%IHN
function [AuO,cor]=scF(Au,FDTot,DTot,IU,W,Tb,Nsamp,ncp,nzc,N,k,Nt,Np,t,TP)
clA=ones(size(Au));
% TP=(N/k*Nsamp+ncp+nzc*Nsamp)*Tb;
for i=1:Nt
    
    for iu=1:Np
        for ij=1:Nt
            if(iu~=IU)
                if(FDTot(iu)-FDTot(IU)<W)
                    if(abs(DTot(iu)-DTot(IU))<TP)
                        if( DTot(iu)>DTot(IU))
                            A=(TP-(DTot(iu)-DTot(IU)))*length(t)/Tb;
                            clA(i,1:floor(A))=zeros(1,floor(A));
                            
                        elseif( DTot(iu)<DTot(IU))
                            Al=clA(i,floor(DTot(IU)*length(1)/Tb):end);
                            clA(i,floor(DTot(IU)*length(1)/Tb):end)=zeros(1,length(Al));
                        end
                    end
                end
            end
        end
    end
end
AuO=Au(1,:);
for i=2:Nt
    Ind1=clA(1,:)==0;
    Ind2=clA(i,:)>0;
    Ind=Ind1.*Ind2; 
    if(sum(Ind)>0)
        AuO(Ind)=Au(i,Ind);
    end
    cor(i)=sum(Ind);
end
cor(1)=length(clA(1,:))-sum(cor(2:end));
cor=cor/sum(clA(1,:));
