
%IHN
%IHN
clc
clear all
close all

Pt=10^((21-30)/10);
alf=1;
Pc=0.010;

InSynT=1;
InSynF=1;
La_Rg=5.3 ;
lambda=La_Rg;
Fm=200;
N=50;%bits/pack
Leng_pack=0.5;
Es=1*Pc;
E0=1000;
Tw=10;
Rc1=50;
Rc2=3000;
Esyt=Pc*3;
Dsyt=3;
Esytf=Pc*5;
Dsytf=5;
FFh= [100:250:3000];
nFFh= [100:100:3000];
Nd=20000;
FasD = sqrt((Rc2^2-Rc1^2)*rand(Nd,1)+Rc1^2);
Ind=zeros(Nd,length(nFFh));
for io=1:length(nFFh)
    Ind(:,io)=  (nFFh(io)-50<FasD & FasD<nFFh(io)+50);
    sIn(io)=sum( Ind(:,io));
end


% kg1=1;
% kg2=7;
% kg=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
H=0;

Nt=1;
co=0;
for lambda=La_Rg
    co=co+1;
    X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H)];
    load([X,'.mat'])
    %     SR111(co)=mean(SucRate);
    NSR111(co,:)=nanmean(nNSucRate);
end

EE=0; 
for j=1
    SR=NSR111 ;
    EE=  La_Rg.*N./((1./SR).*(Leng_pack*(alf*Pt+Pc))+(1./SR-1).*Pc*Tw) ;
end  
  
figure(1)
plot(nFFh ,EE/1000)
hold on

DD=0;
for j=1
    SR=NSR111 ;
    DD=  min(20,(1./SR-1).*(Leng_pack+Tw)+Leng_pack );
end
figure(2)
plot(nFFh ,DD)
hold on
 
LT=0;
for j=1
    SR=NSR111 ;
    LT=  E0./(Es+(1./SR).*Leng_pack*(alf*Pt+Pc)+(1./SR-1).*Pc*Tw);
end

figure(3)
plot(nFFh ,LT)
hold on

kSR=0;
for j=1
    SR=NSR111 ;
    kSR=  SR;
end  

figure(4)
plot(nFFh ,kSR)
hold on

%%%--------------------------

Nt=2;
co=0;
for lambda=La_Rg
    co=co+1;
    X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H)];
    load([X,'.mat'])
    %     SR111(co)=mean(SucRate);
    NSR211(co,:)=nanmean(nNSucRate);
end


EE=0;
for j=1
    SR=NSR211 ;
    EE=  La_Rg.*N./((1./SR).*(2*Leng_pack*(alf*Pt +Pc)+Pc*Leng_pack)+(1./SR-1).*Pc*Tw) ;
end

figure(1)
plot(nFFh ,EE/1000)
hold on


DD=0; 
for j=1
    SR=NSR211 ;
    NSR=NSR111 ;
    Ps=(1-(1-NSR).^2);   
    DD=  min(20,(1./SR-1).*(3*Leng_pack+Tw)+ NSR./Ps.*Leng_pack+(1-NSR).*(NSR)./Ps.*3*Leng_pack);
end
figure(2)
plot(nFFh ,DD)
hold on


LT=0;
for j=1
    SR=NSR211 ;
    LT=  E0./(Es+(1./SR).*(2*Leng_pack*(alf*Pt +Pc)+Pc*Leng_pack)+(1./SR-1).*Pc*Tw);
end

figure(3)
plot(nFFh ,LT)
hold on


kSR=0;
for j=1
    SR=NSR211 ;
    kSR=  SR;
end

figure(4)
plot(nFFh ,kSR)
hold on

%%%--------------------------

Nt=3;
co=0;
for lambda=La_Rg
    co=co+1;
    X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H)];
    load([X,'.mat'])
    %     SR111(co)=mean(SucRate);
    NSR311(co,:)=nanmean(nNSucRate);
end


EE=0; 
for j=1
    SR=NSR311 ;
    EE=  La_Rg.*N./((1./SR).*(3*Leng_pack*(alf*Pt +Pc)+2*Pc*Leng_pack)+(1./SR-1).*Pc*Tw) ;
end 

figure(1)
plot(nFFh ,EE/1000)
hold on

DD=0;
for j=1
    SR=NSR311 ;
    NSR=NSR111 ;
    Ps=(1-(1-NSR).^3);
    DD=  min(20,(1./SR-1).*(5*Leng_pack+Tw)+ NSR./Ps.*Leng_pack+(1-NSR).*(NSR)./Ps.*3*Leng_pack+(1-NSR).^2.*(NSR)./Ps.*5*Leng_pack);
end 
figure(2) 
plot(nFFh ,DD)
hold on



LT=0;
for j=1
    SR=NSR311 ;
    LT=  E0./(Es+(1./SR).*(3*Leng_pack*(alf*Pt +Pc)+2*Pc*Leng_pack)+(1./SR-1).*Pc*Tw);
end

figure(3)
plot(nFFh ,LT)
hold on


kSR=0;
for j=1
    SR=NSR311 ;
    kSR=  SR;
end

figure(4)
plot(nFFh ,kSR)
hold on
%%%--------------------------
H=1;
 
Nt=2;
co=0;
for lambda=La_Rg
    co=co+1;
    X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H)];
    load([X,'.mat'])
    %     SR111(co)=mean(SucRate);
    NSR2111(co,:)=nanmean(nNSucRate);
end


EE=0;
for j=1
    SR=NSR2111 ;
    EE=  La_Rg.*N./((1./SR).*(2*Leng_pack*(alf*Pt/2 +Pc)+Pc*Leng_pack)+(1./SR-1).*Pc*Tw) ;
end

figure(1)
plot(nFFh ,EE/1000)
hold on


DD=0; 
for j=1
    SR=NSR2111 ;
    NSR=NSR111 ;
    Ps=(1-(1-NSR).^2);   
    DD=  min(20,(1./SR-1).*(3*Leng_pack+Tw)+ NSR./Ps.*Leng_pack+(1-NSR).*(NSR)./Ps.*3*Leng_pack);
end
figure(2)
plot(nFFh ,DD)
hold on


LT=0;
for j=1
    SR=NSR2111 ;
    LT=  E0./(Es+(1./SR).*(2*Leng_pack*(alf*Pt/2 +Pc)+Pc*Leng_pack)+(1./SR-1).*Pc*Tw);
end

figure(3)
plot(nFFh ,LT)
hold on


kSR=0;
for j=1
    SR=NSR2111 ;
    kSR=  SR;
end

figure(4)
plot(nFFh ,kSR)
hold on

%%%--------------------------

Nt=3;
co=0;
for lambda=La_Rg
    co=co+1;
    X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H)];
    load([X,'.mat'])
    %     SR111(co)=mean(SucRate);
    NSR3111(co,:)=nanmean(nNSucRate);
end


EE=0;
for j=1
    SR=NSR3111 ;
    EE=  La_Rg.*N./((1./SR).*(3*Leng_pack*(alf*Pt/3 +Pc)+2*Pc*Leng_pack)+(1./SR-1).*Pc*Tw) ;
end

figure(1)
plot(nFFh ,EE/1000)
hold on

DD=0;
for j=1
    SR=NSR3111 ;
    NSR=NSR111 ;
    Ps=(1-(1-NSR).^3);
    DD=  min(20,(1./SR-1).*(5*Leng_pack+Tw)+ NSR./Ps.*Leng_pack+(1-NSR).*(NSR)./Ps.*3*Leng_pack+(1-NSR).^2.*(NSR)./Ps.*5*Leng_pack);
end 
figure(2) 
plot(nFFh ,DD)
hold on



LT=0;
for j=1
    SR=NSR3111 ;
    LT=  E0./(Es+(1./SR).*(3*Leng_pack*(alf*Pt/3 +Pc)+2*Pc*Leng_pack)+(1./SR-1).*Pc*Tw);
end

figure(3)
plot(nFFh ,LT)
hold on


kSR=0;
for j=1
    SR=NSR3111 ;
    kSR=  SR;
end

figure(4)
plot(nFFh ,kSR)
hold on
%%%--------------------------
 
 
%%%--------------------------
H=2;
Nt=3;
co=0;
for lambda=La_Rg
    co=co+1;
    X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF),'H',num2str(H)];
    load([X,'.mat'])
    %     SR111(co)=mean(SucRate);
    NSR4112(co,:)=nanmean(nNSucRate);
end
 

EE=0;
for j=1:length(nFFh)
    SR=NSR4112(j) ;
    if(j<15)
        EE(j)=  La_Rg.*N./((1./SR).*(Leng_pack*(alf*Pt+Pc))+(1./SR-1).*Pc*Tw) ;
    else
        EE(j)=  La_Rg.*N./((1./SR).*(Nt*Leng_pack*(alf*Pt +Pc)+(Nt-1)*Pc*Leng_pack)+(1./SR-1).*Pc*Tw) ;
    end
end

figure(1)
plot(nFFh ,EE/1000)
hold on

DD=0;
for j=1:length(nFFh)
    SR=NSR4112(j) ;
    NSR=NSR111(j) ;
    if(j<15)
        DD(j)=  min(20,(1./SR-1).*(Leng_pack+Tw)+Leng_pack );
        
    else
        Ps=(1-(1-NSR).^3);
        DD(j)=  min(20,(1./SR-1).*(5*Leng_pack+Tw)+ NSR./Ps.*Leng_pack+(1-NSR).*(NSR)./Ps.*3*Leng_pack+(1-NSR).^2.*(NSR)./Ps.*5*Leng_pack);
    end
end 
figure(2)
plot(nFFh ,DD)
hold on



LT=0;
for j=1:length(nFFh)
    SR=NSR4112(j) ;
    if(j<15)
        LT(j)=  E0./(Es+(1./SR).*Leng_pack*(alf*Pt+Pc)+(1./SR-1).*Pc*Tw);
    else
        LT(j)=  E0./(Es+(1./SR).*(3*Leng_pack*(alf*Pt +Pc)+2*Pc*Leng_pack)+(1./SR-1).*Pc*Tw);
    end
end

figure(3)
plot(nFFh ,LT)
hold on 
 

kSR=0;
for j=1:1
    SR=NSR4112 ;
    kSR=  SR;
end

figure(4)
plot(nFFh ,kSR)
hold on
 
%%---------------

figure(1)
legend('N=1; TiAs; FrAs','N=2; TiAs; FrAs','N=3; TiAs; FrAs','ad-N=2; TiAs; FrAs','ad-N=3; TiAs; FrAs','Hyb','ad-Hyb' )
grid on
xlabel('Offered load per {\it W}:  (0.5\timesT_p \times g)')
ylabel('Network Energy Efficiency (Bit/Joule)')

figure(2)
legend('N=1; TiAs; FrAs','N=2; TiAs; FrAs','N=3; TiAs; FrAs','ad-N=2; TiAs; FrAs','ad-N=3; TiAs; FrAs','Hyb','ad-Hyb' )
grid on
xlabel('Offered load per {\it W}:  (0.5\timesT_p \times g)')
ylabel('Packet Delay (Sec)')

figure(3)
legend('N=1; TiAs; FrAs','N=2; TiAs; FrAs','N=3; TiAs; FrAs','ad-N=2; TiAs; FrAs','ad-N=3; TiAs; FrAs','Hyb','ad-Hyb' )
grid on
xlabel('Offered load per {\it W}:  (0.5\timesT_p \times g)')
ylabel('Battery Lifetime (\times reporting period)')

figure(4)
legend('N=1; TiAs; FrAs','N=2; TiAs; FrAs','N=3; TiAs; FrAs','ad-N=2; TiAs; FrAs','ad-N=3; TiAs; FrAs','Hyb','ad-Hyb' )
grid on
xlabel('Offered load per {\it W}:  (0.5\timesT_p \times g)')
ylabel('Success Rate  ')
%%------------
 