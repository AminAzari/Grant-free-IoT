%IHN
function FD=FmF(Np,Fm, InSynF,W)
if(InSynF==1)
    FD=Fm*(2*rand(1,Np)-1);
    
elseif(InSynF==2)
    nF=(2*Fm+W)/W;
    if(nF==2) 
        A=[-1,1];
    elseif(nF==3)
        A=[-2,0,2];
    elseif(nF==4)
        A=[-3,-1,1,3];
    elseif(nF==5) 
        A=[-4,-2,0, 2,4];
    elseif(nF==6)
        A=[-5,-3,-1,1,3,5];
    else
        error('nF odd')
    end
    FD1=1+floor(length(A)*rand(1,Np));
    FD=A(FD1)*W/2;
end