function Hc=cheF(MatPos,Fr,t,ncp,nzc,Tb,Fs)
for it=1:length(Fr)
    dR(it)=FunDrN(nzc,Tb,Fs,-Fr(it));
end  
Hc=ones(size(MatPos,1),1);
for i=1:size(MatPos,1) 
    Hc(i,1)=FunCh(MatPos(i,:),Fr,t,ncp,nzc,Tb,Fs,dR);
end



 