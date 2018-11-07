%IHN
function y=FunDrN(L,Tb,Fs,fd)
seq = transpose(lteZadoffChuSeq(1,L));
%  t=[1/Fs:1/Fs:Tb];
As=[];
for i=1:length(seq)
    AW=(seq(i)*exp(1i*2*pi*(fd)*[(i-1)*Tb+1/Fs:1/Fs:i*Tb]));
    As=[As mean(AW)];
end
Ase=[zeros(1,1*L) As zeros(1,1*L)];

maxd=2*L;
P=zeros(1,maxd);
for d=1:maxd 
    A=power(abs(xcorr(Ase(d:d+L-1),seq)),2);
        P(d)   = A(L);
end 
   
M=P./max(P); 
% figure
% plot(M)
[~,y]=max(M);
y=y-(L+1);
% Ac=abs(xcorr(As,Ase));
% [~,pp]=max(Ac/max(Ac)); 
% y=pp-L; 
% if(y>0)
%     y=y+1;
% elseif(y<0)
%     y=y-1;
% end
%  