%IHN
% function test_th(lambda,Nt,InSynT,InSynF)
lambda=1;
Nt=1;
InSynT=1;
InSynF=1;

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
snr_rg=6;
Eth=0.01;
Pt=21-30;
Rc1=50;
Rc2=1100;
% % % InSynT=1;
% % % InSynF=1;
%N/k*Nsamp+ncp+nzc*Nsamp  254
W=2/Tb;
Fm=200;
% % % lambda=1.5;  %%Nt=2,Eth=0.1,lamb=0.7-->rel=0.9
% % % Nt=3;
r1=2;r2=2*Nt;
Leng_pack=(N/k*Nsamp+ncp+nzc*Nsamp)*Tb;
off_load=lambda*Leng_pack;
% off_load_eff=off_load/((2*Fm+W)/W);
% lambda_eff=lambda/((2*Fm+W)/W);
% lambda_eff_pack=lambda/((2*Fm+W)/W)*Nt;
% off_load_eff_pack=off_load/((2*Fm+W)/W)*Nt;

Tperiod=50;%sec
meanNp=1;
disp(['Np=',num2str(meanNp),'%%% Nt=',num2str(Nt),'%%% off_load_pack=',num2str(off_load*Nt),'%%% off_load_eff=',num2str(off_load)])

t=[1/Fs:1/Fs:Tb];

cons=.1;
ia=0;
Nrep=1000;
SucRate=zeros(1,Nrep);
NSucRate=zeros(Nrep,2);
ConAl=zeros(3,Nrep);
for IA=1:Nrep
    ia=ia+1;
    if(mod(ia,50)==0)
    disp(['ItCont/Nrep=',num2str(ia/Nrep)])
    end
    Np=1;
    if(Np==0)
        Np=1;
    end
    D=TiOfFun(Np,t,Tb,r1,r2,Nt,Leng_pack,Tperiod,InSynT);
    FD=FmF(Np,Fm, InSynF,W);
    
    
    %%-TRansmitter
    
    [y_pam_temp,ip,y0,cy0,fas]=FunTxP(N,k,Nsamp, ncp,nzc,fc,Fs,Tb,FD,D,cons,Pt,Rc1,Rc2,W,Nt);
    t=[1/Fs:1/Fs:Tb];
    
    %%%-----Channel
    %embeded in Tx
    
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
        if(Nt>1)
            In=min(nE)<Eth;
        else
            In= nE<Eth;
        end
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
    %%display('%%%%%%%---egc---%%%%%%%%%%%%%%%%%')
    if(sum(Sucn)~=Np)
        com='egc';
        Wcon=1;
        while(Wcon)
            nE=zeros(Nt,Np);
            for pco=1:Np
                if(Suc(pco)==0)
                    nET=decF(FD(pco),D(:,pco),y_pam_temp,ip(pco,:),Fs,Tb,N,k,Nsamp,ncp,nzc,fc,com,FD,D,y0,cy0,pco,W,Nt,Np);
                    nE(:,pco)=nET*ones(Nt,1);
                end
            end
            if(Nt>1)
                In=min(nE)<Eth;
            else
                In= nE<Eth;
            end
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
                X = ['egc,success rate= ',num2str(sum(Sucn)/Np)];
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
            if(Nt>1)
                In=min(nE)<Eth;
            else
                In= nE<Eth;
            end
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
    
   
    FFh= [100:250:3000];
    for io=1:length(FFh)
        Ind=  (FFh(io)-50<fas & fas<FFh(io)+50);
        %disp(['cont=',num2str(sum(ia)),'%%% Sum Ssuc=',num2str(sum(Sucn)/Np)])
        if(sum(Ind)>0)
            NSucRate(IA,io)=sum(Sucn(Ind))/sum(Ind);
        else
            NSucRate(IA,io)=nan;
        end
    end
    SucRate(IA)=sum(Sucn)/Np;
end
nanmean(NSucRate(:,5))
 