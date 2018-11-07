%IHN
clc
clear all
close all
nzc=23;
ncp=2;
N = 50;  %number of bits or symbols
M = 4;    %QPSK constellation
k = log2(M);  %number of bits per symbol
Tb=10*1e-3;%1e-2; %100ps
bit_rate=1/Tb; %100bps
fc=800*1e6;%100e9; %carrier=100GHz
Fs=40/Tb; %sampling rate This is the highest freq that you can see on PSD, 2pi*Fs/2pi=400
%note: signal BW in baseband:1/(Tb*nsamp) =200; in passband:2/Tb=400
%note: here we expect to track CFOs up to 2500 wit 5/Tb=Fs
Ts=k*Tb;
Nsamp=1; 
snr_rg=1;
Eth=0.01;
InSynT=1;
InSynF=1;
%N/k*Nsamp+ncp+nzc*Nsamp  254
W=2/Tb;
Fm=100;
lambda=1.5;  %%Nt=2,Eth=0.1,lamb=0.7-->rel=0.9
Nt=3;
r1=2;r2=2*Nt;  
Leng_pack=(N/k*Nsamp+ncp+nzc*Nsamp)*Tb;
off_load=lambda*Leng_pack;
off_load_eff=off_load/((2*Fm+W)/W);
lambda_eff=lambda/((2*Fm+W)/W);
lambda_eff_pack=lambda/((2*Fm+W)/W)*Nt;
off_load_eff_pack=off_load/((2*Fm+W)/W)*Nt;

Tperiod=50;%sec
meanNp=floor(Tperiod*lambda);
disp(['Np=',num2str(meanNp),'%%% Nt=',num2str(Nt),'%%% off_load_pack=',num2str(off_load*Nt),'%%% off_load_eff=',num2str(off_load)])

t=[1/Fs:1/Fs:Tb];

cons=.1;
ia=0;
Nrep=100;
SucRate=zeros(1,Nrep);
ConAl=zeros(3,Nrep); 
for IA=1:Nrep
    ia=ia+1; 
    Np=poissrnd(meanNp,1,1);
     D=TiOfFun(Np,t,Tb,r1,r2,Nt,Leng_pack,Tperiod,InSynT);
    FD=FmF(Np,Fm, InSynF,W);
    
    %%-TRansmitter
    [y,ip,y0,cy0]=FunTxP(N,k,Nsamp, ncp,nzc,fc,Fs,Tb,FD,D,cons);
    t=[1/Fs:1/Fs:Tb];
    %%-Channel
    [yp,~]=FunTxP(N,k,Nsamp, ncp,nzc,fc,Fs,Tb,0,0,cons);
    n=1/sqrt(2)*[randn(1,length(y)) + 1i*randn(1,length(y))];
    py=sum(yp*yp')/length(yp);
    pn=sum(n*n')/length(n);
    sn=snr_rg-10*log10(py/pn) ;
    y_pam_temp = y +10^(-(sn)/20)*n;  %additive white gaussian noise%+awgn(y,snr_rg(ii));
    
    Suc=zeros(1,Np);
    
    %%display('%%%%--single---%%%%%%%%%%%')
    wc=0;
    Wcon=1;
    while(Wcon)
        nE=zeros(Nt,Np);
        for pco=1:Np
            if(Suc(pco)==0)
                for tco=1:Nt
                    nE(tco,pco)=decF(FD(pco),D(tco,pco),y_pam_temp,ip(pco,:),Fs,Tb,N,k,Nsamp,ncp,nzc,fc);
                end
            end
        end
        In=min(nE)<Eth;
        Sucn=Suc;
        Sucn(In)=ones(1,sum(In));
        dSuc=Sucn-Suc;
        y_pam_temp=remF(y_pam_temp,dSuc,D,FD,Fs,Tb,ip,N,nzc,ncp,Nsamp,Np,fc,k);
        
        Suc=Sucn;
        %disp(['sum dsuc= ',num2str(sum(dSuc))])
        if(sum(dSuc)>0)
            wc=wc+1;
            %disp(['contr= ',num2str(wc)])
%             nEr(:,:,wc)=nE;
        end
        if(  sum(Sucn)==Np)
            Wcon=0;
            %%display('success 100%,single')
        end
        if( sum(dSuc)==0)
            Wcon=0;
            %disp(['single,success rate= ',num2str(sum(Sucn)/Np)])
        end
        %%display('--------------------------------')
    end
    ConAl(1,ia)=sum(Sucn)/Np;
    %%display('%%%%%%%---MRC---%%%%%%%%%%%%%%%%%')
    if(sum(Sucn)~=Np)
        com='mrc';
        Wcon=1;
        while(Wcon)
            nE=zeros(Nt,Np);
            for pco=1:Np
                if(Suc(pco)==0)
                    nET=decF(FD(pco),D(:,pco),y_pam_temp,ip(pco,:),Fs,Tb,N,k,Nsamp,ncp,nzc,fc,com,FD,D,y0,cy0,pco,W,Nt,Np);
                    nE(:,pco)=nET*ones(Nt,1);
                end
            end
            In=min(nE)<Eth;
            Sucn=Suc;
            Sucn(In)=ones(1,sum(In));
            dSuc=Sucn-Suc;
            y_pam_temp=remF(y_pam_temp,dSuc,D,FD,Fs,Tb,ip,N,nzc,ncp,Nsamp,Np,fc,k);
            
            Suc=Sucn;
            %disp(['sum dsuc= ',num2str(sum(dSuc))])
            if(sum(dSuc)>0)
                wc=wc+1;
                %disp(['contr= ',num2str(wc)])
%                 nEr(:,:,wc)=nE;
            end
            
            if(  sum(Sucn)==Np)
                Wcon=0;
                %%display('success 100%,mrc')
            end
            if( sum(dSuc)==0)
                Wcon=0;
                X = ['mrc,success rate= ',num2str(sum(Sucn)/Np)];
                %disp(X)
            end
        end
        %%display('--------------------------------')
        ConAl(2,ia)=sum(Sucn)/Np-ConAl(1,ia);
    else
        ConAl(2,ia)=nan;
    end
    
    
    %%display('%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    
    
    
    %%display('%%%%%%%---SC---%%%%%%%%%%%%%%%%%')
    if(sum(Sucn)~=Np)
        com='sc';
        Wcon=1;
        while(Wcon)
            nE=zeros(Nt,Np);
            for pco=1:Np
                if(Suc(pco)==0)
                    nET=decF(FD(pco),D(:,pco),y_pam_temp,ip(pco,:),Fs,Tb,N,k,Nsamp,ncp,nzc,fc,com,FD,D,y0,cy0,pco,W,Nt,Np);
                    nE(:,pco)=nET*ones(Nt,1);
                end
            end
            In=min(nE)<Eth;
            Sucn=Suc;
            Sucn(In)=ones(1,sum(In));
            dSuc=Sucn-Suc;
            y_pam_temp=remF(y_pam_temp,dSuc,D,FD,Fs,Tb,ip,N,nzc,ncp,Nsamp,Np,fc,k);
            
            Suc=Sucn;
            %disp(['sum dsuc= ',num2str(sum(dSuc))])
            if(sum(dSuc)>0)
                %disp(['contr= ',num2str(wc)])
%                 nEr(:,:,wc)=nE;
            end
            
            if(  sum(Sucn)==Np)
                Wcon=0;
                %%display('success 100%,sc')
            end
            if( sum(dSuc)==0)
                Wcon=0;
                %disp(['sc,success rate= ',num2str(sum(Sucn)/Np)])
            end
        end
        %%display('--------------------------------')
        ConAl(3,ia)=sum(Sucn)/Np-ConAl(2,ia)-ConAl(1,ia);
    else
        ConAl(3,ia)=nan;
    end
    %%display('%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    SucRate(IA)=sum(Sucn)/Np;
    disp(['cont=',num2str(sum(ia)),'%%% Sum Ssuc=',num2str(sum(Sucn)/Np)])
end
mean(SucRate)

