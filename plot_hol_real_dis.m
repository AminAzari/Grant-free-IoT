
%IHN
%IHN
clc
clear all
close all

Pt=10^((21-30)/10);
alf=1/0.5;
Pc=0.05;
InSynT=1;
InSynF=1;
La_Rg=0.1:.4:1.7;
lambda=La_Rg;
Fm=200;
N=50;%bits/pack
Leng_pack=0.5;
Es=1*Pc;
E0=1000;
Tw=2;
Esyt=Pc*3;
Dsyt=3;
Esytf=Pc*5;
Dsytf=5;
kg=9;

Nt=1;
co=0;
for lambda=La_Rg
    co=co+1;
    X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF)];
    load([X,'.mat'])
    SR111(co)=mean(SucRate);NSR111(co,:)=nanmean(NSucRate);
    if(co==1)
        Aa=[];
        for i=1:length(La_Rg)
            Aa=[Aa;nanmean(NSucRate)];
        end
        SRn(1:length(La_Rg))=mean(SucRate)*ones(1,length(La_Rg));NSRn(1:length(La_Rg),:)=Aa;
    end
end

SR111=mean(NSR111(:,kg),2)';
SR=SR111;

figure(1) 
EE=La_Rg.*N./((1./SR).*(Leng_pack*(alf*Pt+Pc)+(1./SR-1).*Pc*Tw));
plot(La_Rg ,EE/1000)
hold on
 
figure(2)
DD=((1./SR-1).*(Leng_pack+Tw)+Leng_pack );
plot(La_Rg ,DD)
hold on

figure(3)
LT=E0./(Es+(1./SR).*Leng_pack*(alf*Pt+Pc)+(1./SR-1).*Pc*Tw);
plot(La_Rg ,LT)
hold on

figure(4)
plot(La_Rg ,SR)
hold on
%%%--------------------------

Nt=2;
co=0;
for lambda=La_Rg
    co=co+1;
    X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF)];
    load([X,'.mat'])
    SR211(co)=mean(SucRate);NSR211(co,:)=nanmean(NSucRate);
end
SR211=NSR211(:,kg)';
SR=SR211;
figure(1)
EE=La_Rg.*N./((1./SR).*(Nt*Leng_pack*(alf*Pt/Nt+Pc)+(1./SR-1).*Pc*Tw));
plot(La_Rg ,EE/1000)
hold on

 Ps=2.*SR111-SR111.^2;

figure(2)
DD=((1./SR-1).*(3*Leng_pack+Tw)+ SR111./Ps.*Leng_pack+(1-SR111).*(SR111)./Ps.*3*Leng_pack);
plot(La_Rg ,DD) 
hold on

figure(3)
LT=E0./(Es+(1./SR).*Nt*Leng_pack*(alf*Pt/Nt+Pc)+(1./SR-1).*Pc*Tw);
plot(La_Rg ,LT)
hold on

figure(4)
plot(La_Rg ,SR)
hold on
%%%--------------------------

% Nt=3;
% co=0;
% for lambda=La_Rg
%     co=co+1;
%     X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF)];
%     load([X,'.mat'])
%     SR311(co)=mean(SucRate);NSR311(co,:)=nanmean(NSucRate);
% end
% 
% SR311=NSR311(:,kg)';
% SR=SR311;
% figure(1)
% EE=La_Rg.*N./((1./SR).*(Nt*Leng_pack*(alf*Pt/Nt+Pc)+(1./SR-1).*Pc*Tw));
% plot(La_Rg ,EE/1000)
% hold on
% 
%   Ps=1-(1-SR111).^3; 
% 
% figure(2)
% DD=((1./SR-1).*(5*Leng_pack+Tw)+ SR111./Ps.*Leng_pack+(1-SR111).*(SR111)./Ps.*3*Leng_pack+(1-SR111).^2.*(SR111)./Ps.*5*Leng_pack);
% plot(La_Rg ,DD)
% hold on
% 
% figure(3)
% LT=E0./(Es+(1./SR).*Nt*Leng_pack*(alf*Pt/Nt+Pc)+(1./SR-1).*Pc*Tw);
% plot(La_Rg ,LT)
% hold on
% 
% figure(4)
% plot(La_Rg ,SR)
% hold on
%%%--------------------------

InSynT=2;
InSynF=1;
Nt=1;
co=0;
for lambda=La_Rg
    co=co+1;
    X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF)];
    load([X,'.mat'])
    SR121(co)=mean(SucRate);NSR121(co,:)=nanmean(NSucRate);
end

SR121=NSR121(:,kg)';
SR=SR121;
figure(1)
EE=La_Rg.*N./(Esyt+(1./SR).*(Leng_pack*(alf*Pt+Pc)+(1./SR-1).*Pc*Tw));
plot(La_Rg ,EE/1000)
hold on

figure(2)
DD=(Dsyt+(1./SR-1).*(Leng_pack+Tw)+Leng_pack );
plot(La_Rg ,DD)
hold on

figure(3)
LT=E0./(Esyt+Es+(1./SR).*Leng_pack*(alf*Pt+Pc)+(1./SR-1).*Pc*Tw);
plot(La_Rg ,LT)
hold on

figure(4)
plot(La_Rg ,SR)
hold on
%%%--------------------------

% InSynT=2;
% InSynF=1;
% Nt=2;
% co=0;
% for lambda=La_Rg
%     co=co+1;
%     X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF)];
%     load([X,'.mat'])
%     SR221(co)=mean(SucRate);NSR221(co,:)=nanmean(NSucRate);
% end
% 
% SR221=NSR221(:,kg)';
% SR=SR221;
% figure(1)
% EE=La_Rg.*N./(Esyt+(1./SR).*(Nt*Leng_pack*(alf*Pt/Nt+Pc)+(1./SR-1).*Pc*Tw));
% plot(La_Rg ,EE/1000)
% hold on
% 
% Ps=1-(1-SR111).^2;
% figure(2)
% DD=(Dsyt+(1./SR-1).*(3*Leng_pack+Tw)+ SR121./Ps.*Leng_pack+(1-SR121).*(SR121)./Ps.*3*Leng_pack);
% plot(La_Rg ,DD)
% hold on
% 
% figure(3)
% LT=E0./(Esyt+Es+(1./SR).*Nt*Leng_pack*(alf*Pt/Nt+Pc)+(1./SR-1).*Pc*Tw);
% plot(La_Rg ,LT)
% hold on
% 
% figure(4)
% plot(La_Rg ,SR)
% hold on
% %%%--------------------------
InSynT=2;
InSynF=2;
Nt=1;
co=0;
for lambda=La_Rg
    co=co+1;
    X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF)];
    load([X,'.mat'])
    SR122(co)=mean(SucRate);NSR122(co,:)=nanmean(NSucRate);
end

SR122=NSR122(:,kg)';
SR=SR122;
figure(1)
EE=La_Rg.*N./(Esytf+(1./SR).*(Leng_pack*(alf*Pt+Pc)+(1./SR-1).*Pc*Tw));
plot(La_Rg ,EE/1000)
hold on

figure(2)
DD=(Dsytf+(1./SR-1).*(Leng_pack+Tw)+Leng_pack );
plot(La_Rg ,DD)
hold on

figure(3)
LT=E0./(Esytf+Es+(1./SR).*Leng_pack*(alf*Pt+Pc)+(1./SR-1).*Pc*Tw);
plot(La_Rg ,LT)
hold on

figure(4)
plot(La_Rg ,SR)
hold on
%%%--------------------------

 
InSynT=2;
InSynF=2;
Nt=2;
co=0;
for lambda=La_Rg
    co=co+1;
    X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF)];
    load([X,'.mat'])
    SR222(co)=mean(SucRate);NSR222(co,:)=nanmean(NSucRate);
end

SR222=NSR222(:,kg)';
SR=SR222;
figure(1)
EE=La_Rg.*N./(Esytf+(1./SR).*(Nt*Leng_pack*(alf*Pt/Nt+Pc)+(1./SR-1).*Pc*Tw));
plot(La_Rg ,EE/1000)
hold on

Ps=1-(1-SR111).^2;
figure(2)
DD=(Dsytf+(1./SR-1).*(3*Leng_pack+Tw)+ SR122./Ps.*Leng_pack+(1-SR122).*(SR122)./Ps.*3*Leng_pack);
plot(La_Rg ,DD)
hold on  

figure(3)
LT=E0./(Esytf+Es+(1./SR).*Nt*Leng_pack*(alf*Pt/Nt+Pc)+(1./SR-1).*Pc*Tw);
plot(La_Rg ,LT)
hold on

figure(4)
plot(La_Rg ,SR)
hold on
%%%--------------------------
Nres=10;
Tpr=2;
InSynT=2;
InSynF=2;
Nt=1;
co=0;
for lambda=La_Rg
    co=co+1;
    Pss(co)=FunPs(Nres,lambda,Tpr);
    Del=Tpr/2+max(Tpr,1/(11/Tpr-lambda));
end

figure(1)
EE=La_Rg.*N./(Esytf+Pc*Del+(1./SRn).*(1./Pss).*(Leng_pack/5*(alf*Pt+Pc)+Pc*1)+(1./SRn)*(Leng_pack*(alf*Pt+Pc)+2*Pc*Tpr));
plot(La_Rg ,EE/1000)
hold on

figure(2)
DD=(Dsytf+Del +(1./Pss-1).*Tpr++(1./SRn)*(Leng_pack+Pc*2*Tpr) );
plot(La_Rg ,DD)
hold on

figure(3)
LT=E0./(Esytf+Pc*Del+(1./SRn).*(1./Pss).*(Leng_pack/5*(alf*Pt+Pc)+Pc*1)+(1./SRn)*(Leng_pack*(alf*Pt+Pc)+Pc*2*Tpr));
plot(La_Rg ,LT)
hold on
%%---------------

figure(1)
legend('N=1; TiAs; FrAs','N=2; TiAs; FrAs','N=3; TiAs; FrAs','N=1; TiSy; FrAs','N=2; TiSy; FrAs','N=1; TiSy; FrSy','N=2; TiSy; FrSy','Granted (Synch)')
grid on
xlabel('Arrival rate of packets (per second)')
ylabel('Network Energy Efficiency (Bit/Joule)')

figure(2)
legend('N=1; TiAs; FrAs','N=2; TiAs; FrAs','N=3; TiAs; FrAs','N=1; TiSy; FrAs','N=2; TiSy; FrAs','N=1; TiSy; FrSy','N=2; TiSy; FrSy','Granted (Synch)')
grid on
xlabel('Arrival rate of packets (per second)')
ylabel('Packet Delay (Sec)')

figure(3)
legend('N=1; TiAs; FrAs','N=2; TiAs; FrAs','N=3; TiAs; FrAs','N=1; TiSy; FrAs','N=2; TiSy; FrAs','N=1; TiSy; FrSy','N=2; TiSy; FrSy','Granted (Synch)')
grid on
xlabel('Arrival rate of packets (per second)')
ylabel('Battery Lifetime (\times reporting period)')

figure(4)
legend('N=1; TiAs; FrAs','N=2; TiAs; FrAs','N=3; TiAs; FrAs','N=1; TiSy; FrAs','N=2; TiSy; FrAs','N=1; TiSy; FrSy','N=2; TiSy; FrSy' )
grid on
xlabel('Arrival rate of packets (per second)')
ylabel('Success Rate  ')
%%------------
 

