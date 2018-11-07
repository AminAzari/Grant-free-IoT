%%--------------------------timing recovery
function [co,chc,chcr,coo,coE]=FunTim_dm_EFilAl(y_dqpsk_temp,t,fc,Tb,Fs,ncp,nzc,N,k, Nsamp,Fr,tune,a1,a2,a3,th)
y_rt=[];
j=0;

for i=1:length(t):length(y_dqpsk_temp)
    j=j+1;
    AAk0=y_dqpsk_temp(i:i+length(t)-1).*exp(-1i*2*pi*(fc)*[(j-1)*Tb+1/Fs:1/Fs:j*Tb]);
    y_rt(1,(j-1)*length(t)+1:j*length(t))=(AAk0) ; %1+floor(20*rand(1,20)
    for w=1:length(Fr)
        AAk=y_dqpsk_temp(i:i+length(t)-1).*exp(-1i*2*pi*(fc+Fr(w))*[(j-1)*Tb+1/Fs:1/Fs:j*Tb]);
        y_rt(w+1,(j-1)*length(t)+1:j*length(t))= (AAk) ;%(1+floor(20*rand(1,20)))
    end
end



maxd=length(y_rt(1,:))-(ncp+nzc*Nsamp+N/k*Nsamp)*length(t)-5;

P=zeros(length(Fr)+1,maxd);
P0=zeros(length(Fr)+1,maxd);
received_frame=y_rt;

Ase= (transpose(lteZadoffChuSeq(1,nzc)));
Ase= rectpulse(Ase,Nsamp);
Ase=rectpulse(Ase,length(t));

for d=1:maxd
    P(1,d)   = power(abs(sum(conj(received_frame(1,d:d+ncp*length(t)-1)).*(received_frame(1,d+(N/k*Nsamp+nzc*Nsamp)*length(t):d+(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)-1)))),a3);
    for u=1:length(Fr)
        P(u+1,d)   = power(abs(sum(conj(received_frame(u+1,d+ncp*length(t):d+(nzc*Nsamp+ncp)*length(t)-1)).*Ase)),a1);
        P0(u+1,d)   = power(abs(sum(conj(received_frame(u+1,d:d+ncp*length(t)-1)).*(received_frame(u+1,d+(N/k*Nsamp+nzc*Nsamp)*length(t):d+(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)-1)))),a2);
    end
end

for d=1:maxd
    for u=1:length(Fr)
        P(u+1,d)   = P(u+1,d)*P0(u+1,d)*P(1,d);
    end
end

% % close all
% dp=1:1:maxd;
% figure()
% MM1=P(2,:)/max(P(2,:));
% MM2=P(3,:)/max(P(3,:));
% MM3=P(4,:)/max(P(4,:));
% MM4=P(5,:)/max(P(5,:));
% MM5=P(6,:)/max(P(6,:));
% subplot(4,1,1);plot(dp,MM1);
% % xlim([1 1600])
% subplot(4,1,2);plot(dp,MM2);
% % xlim([1 1600])
% subplot(4,1,3);plot(dp,MM3);
% % % xlim([1 1600])
% subplot(4,1,4);plot(dp,MM4);
% % xlim([1 1600])
% subplot(5,1,5);plot(dp,MM5);
% % xlim([1 1600])
% grid on;


co=zeros(length(Fr),1);
for u=1:length(Fr)
    M1(u+1,:)=P(u+1,:)/max(P(u+1,:));
end
MaxM=zeros(length(Fr),tune);
for u=1:length(Fr)
    M=M1(u+1,:);
    sM=sort(M);
    [V_max_M,P_max_M]=findpeaks(M,'MinPeakHeight',sM(floor(th*length(sM))));
    if length(V_max_M)>tune
        [~,so_V_max_M_Info]=sort(V_max_M);
        so_P_max_M=P_max_M(so_V_max_M_Info);
        P_max_M=so_P_max_M(end-tune+1:end);
        %          V_max_M=so_V_max_M(end-length(Fr)+1:end);
    end
    MaxM(u,1:length(P_max_M))=P_max_M;
end 
%
co=ComF(MaxM);
chc=cheF(co,Fr,t,ncp,nzc,Tb,Fs);
rankk=zeros(size(co,1),length(Fr));
for iy=1:size(co,1)
    for iu=1:length(Fr)
        M=M1(iu+1,:);
        rankk(iy,iu)=FunRaJ(co(iy,iu),M,1);%floor(length(t)-1)
    end
end
chcr=10*log10(prod(rankk,2));


chd=chc==0;
A=chd.*chcr;
[sA,sAI]=sort(A,'descend');
sf=find(sA==0);

if(size(co,1)>1 && sum(sf)>0)
    coo=co(sAI(1:sf(1)-1),:);
    coE=sA(1:sf(1)-1);
else
    coo=co;
    coE=sA;
end



