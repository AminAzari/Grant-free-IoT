%IHN
%Simple Matlab/Octave code for single polarity PAM
clc
clear
close all
clear all
nzc=49;
ncp=5;
N = 4000;  %number of bits or symbols
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
snr_rg=50;
%N/k*Nsamp+ncp+nzc*Nsamp  254
FD=[1*fd  2*fd 2.3*fd];
D=[ 1     1      3000
    3000  7000   7000
    10000 15000  15000];
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
    ia=ia+1
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
        
        %%----decoding
        %spooM=PerPo;
        spooM=D;
        sFr=FD;
        cfo=sFr(1);
        yk(1,:)=y_pam_temp(D(1,1)*length(t)+1:(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)+D(1,1)*length(t));
        yk(2,:)=y_pam_temp(D(2,1)*length(t)+1:(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)+D(2,1)*length(t));
        ak=y_pam_temp(D(3,1)*length(t)+1:(N/k*Nsamp+ncp+nzc*Nsamp)*length(t)+D(3,1)*length(t));
        yyk(1,:)= FuncDecP_inf(t,yk(1,:),fc,Fs,Tb,cfo,ncp,nzc,Nsamp);
        yyk(2,:)= FuncDecP_inf(t,yk(2,:),fc,Fs,Tb,cfo,ncp,nzc,Nsamp);
        aak= FuncDecP_inf(t,ak,fc,Fs,Tb,cfo,ncp,nzc,Nsamp);
        
        rak=real(aak);
        Nl=length(aak);
        sy=sum(rak.^2)/Nl;
        if1=real(aak)-real(yyk(1,:)); 
        if2=real(aak)-real(yyk(2,:));
        sr1=sum(if1.^2)/Nl;
        sr2=sum(if2.^2)/Nl; 
        a1=sy/(sy*(1+sr1/sr2)+sr1);
        a2=sy/(sy*(1+sr2/sr1)+sr2);
        yyk(3,:)=(a1*yyk(1,:)+a2*yyk(2,:))/(a1+a2);
        if3=real(aak)-real(yyk(3,:));
        sr3=sum(if3.^2)/Nl;
        aA(ia,1:2)=[a1 a2]; 
        [sy/sr1 sy/sr2 sy/sr3]
          
        cv=3; 
        RR=[2:0.1:3];
        for iy=RR
            cv=cv+1;
            yk(cv,:)=(1*yk(1,:)+iy*yk(2,:))/(1+iy);
        end
        
         
        
        for iij=1:2+length(RR)+1;
            y_pam_temp_iij=yk(iij,:);
            estBit_noncoh= FuncDecP(t,y_pam_temp_iij,fc,Fs,Tb,cfo,ncp,nzc,Nsamp);
            %--counting errorsm
            ipb=ip(1,:)>mean(ip(1,:));
            nErr_pam_noncoh(ia,iij,ii) = size(find([ipb - estBit_noncoh]),2);  %counting the number of errors
        end
        
    end
end
mean(nErr_pam_noncoh/N)
plot(10*log10(mean(nErr_pam_noncoh/N)))
mean(aA(:,2)./aA(:,1))
[a1,a2]=min(mean(nErr_pam_noncoh/N));
RR(a2-2)
