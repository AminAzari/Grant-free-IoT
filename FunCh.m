function y=FunCh(PoS,Fr,t,ncp,nzc,Tb,Fs,dR)
y=0; 
for u=1:length(Fr)
    po=PoS(u);
    for ij=1:length(Fr)
        if(ij~=u)
            poij=PoS(ij);
            mpo=poij-po;
            ffd=length(t)*(dR(ij) - dR(u) );
            if(abs(mpo-ffd)<length(t)/3)
                y=y+1;
            end
        end
    end
    
end
