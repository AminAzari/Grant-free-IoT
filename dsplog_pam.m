%IHN
%Simple Matlab/Octave code for single polarity PAM
clc
clear
close all
clear all
nzc=49;
ncp=5;
N = 400;  %number of bits or symbols
M = 4;  %QPSK constellation
k = log2(M);  %number of bits per symbol
tune=4;
Tb=10*1e-3;%1e-2; %100ps
bit_rate=1/Tb; %100bps
fc=800*1e6;%100e9; %carrier=100GHz
Fs=40/Tb; %sampling rate This is the highest freq that you can see on PSD, 2pi*Fs/2pi=400
%note: signal BW in baseband:1/(Tb*nsamp) =200; in passband:2/Tb=400
%note: here we expect to track CFOs up to 2500 wit 5/Tb=Fs
Ts=k*Tb;
fd=50;

Nsamp=1;
dx=50;
snr_rg=-2;
%N/k*Nsamp+ncp+nzc*Nsamp  254
FD=[1*fd 2*fd  ];
D=[ 2 2 ];
%ne=zeros(1,length(snr_rg));
%neEx=zeros(1,length(snr_rg));
%np=zeros(1,length(snr_rg));

a1=30;
a2=0;
a3=0;
ia=0;
th=0.99;
cons=.1;

Nrep=10;
for IA=1:Nrep
    ia=ia+1;
    for ii = 1:length(snr_rg)
        %ii
        %%-TRansmitter
        [y,ip]=FunTxP(N,k,Nsamp, ncp,nzc,fc,Fs,Tb,FD,D,dx,cons);
        t=[1/Fs:1/Fs:Tb];
        %%-Channel
        [yp,~]=FunTxP(5*N,k,Nsamp, ncp,nzc,fc,Fs,Tb,0,0,0,cons);
        n=1/sqrt(2)*[randn(1,length(y)) + 1i*randn(1,length(y))];
        py=sum(yp*yp')/length(yp);
        pn=sum(n*n')/length(n);
        sn=snr_rg(ii)-10*log10(py/pn) ;
        y_pam_temp = y +10^(-(sn)/20)*n;  %additive white gaussian noise%+awgn(y,snr_rg(ii));
        
        %%Receiver 
        %load data.mat 
        %%---finding CFOs
        %periodogram(y_pam_temp,[],Fs,'centered')
        [pk,p1]=periodogram(y_pam_temp,[],Fs,'centered');
        HX=p1/pi;
        HY=10*log10(pk);
        mHY=max(HY); 
        [~, po]=findpeaks(HY,'MinPeakDistance',1,'MinPeakHeight',mHY-5);
        Fr=(HX(po)*Fs/2);
        
        clear ij
        %load data.mat
        [pos,chc,chcj,Coo,CooC]=FunTim_dm_EFilAl([   y_pam_temp ],t,fc,Tb,Fs,ncp,nzc,N,k,Nsamp,Fr,tune,a1,a2,a3,th );%all javab
        poo=(pos-1)/length(t);
        spoo=sort(poo);
        [~,sInfo]=sort(poo(1,:));
        aw=zeros(1,size(spoo,1));
        for iq=1:size(spoo,1)
            aw(iq)=size(find([spoo(iq,:) - D]),2);
        end
        [ne(ia,ii),~]=min(aw);
        
        coo=(Coo-1)/length(t);
        SCoo=coo(:,sInfo);
        aw=zeros(1,size(SCoo,1));
        for iq=1:size(SCoo,1)
            aw(iq)=size(find([SCoo(iq,:) - D]),2);
        end
        
        [ne2(ia,ii),mak]=min(aw);
        PerP=Coo(mak,:);
        PerPo=(PerP-1)/length(t);
        pooM=Coo(1,:);
        posM=(pooM-1)/length(t);
        spooM=posM(sInfo);
        aw0=size(find([spooM - D]),2);
        ne0(ia,ii)=min(aw0);
        sFr=Fr(sInfo);
        
        %yk(1,:)=y_pam_temp(D(1,1)*length(t)+1:(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)+D(1,1)*length(t));
        %yk(3,:)=y_pam_temp(D(2,1)*length(t)+1:(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)+D(2,1)*length(t));
        %yk(2,:)=y_pam_temp(D(1,2)*length(t)+1:(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)+D(1,2)*length(t));
        %yk(4,:)=y_pam_temp(D(2,2)*length(t)+1:(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)+D(2,2)*length(t));
        
        
        %%----decoding
        %spooM=PerPo;
        spooM=D;
        sFr=FD;
        for iij=1:length(Fr)
            y_pam_temp_iij=y_pam_temp(spooM(iij)*length(t)+1:(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)+spooM(iij)*length(t));
            %y_pam_temp_iij=yk(iij,:)+yk(iij+2,:);
            cfo=sFr(iij);
            estBit_noncoh= FuncDecP(t,y_pam_temp_iij,fc,Fs,Tb,cfo,ncp,nzc,Nsamp);
            %--counting errors
            ipb=ip(iij,:)>mean(ip(iij,:));
            nErr_pam_noncoh(ia,iij,ii) = size(find([ipb - estBit_noncoh]),2);  %counting the number of errors
%             if(nErr_pam_noncoh(ia,iij,ii)>0)
%                 [sFr(iij) FD(iij) spooM(iij) D(iij)]
%             end
        end
        
    end
end
[mean(ne)/length(Fr),mean(ne0)/length(Fr),mean(ne2)/length(Fr)]
mean(nErr_pam_noncoh/N)

