
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
La_Rg=1.3;
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
kg=7;
% lambda=La_Rg;

co=0;
for Nt=1:5
    co=co+1;
    X=['La',num2str(lambda),'Nt',num2str(Nt),'InSynT',num2str(InSynT),'InSynF',num2str(InSynF)];
    load([X,'.mat'])
    SR111(co)=mean(SucRate);NSR111(co,:)=nanmean(NSucRate);
end

SR111=mean(NSR111(:,kg),2)';
SR=SR111;
Nt_Rg=[1:5];
EE=La_Rg.*N./((1./SR).*(Nt_Rg.*Leng_pack.*(alf.*Pt./Nt_Rg+Pc)+(1./SR-1).*Pc*Tw));
plot(Nt_Rg ,EE/1000)
hold on
 
figure(2)
LT=E0./(Es+(1./SR).*Nt_Rg.*Leng_pack.*(alf*Pt./Nt_Rg+Pc)+(1./SR-1).*Pc*Tw);
plot(Nt_Rg ,LT)
hold on

figure(3)
plot(Nt_Rg ,SR)
hold on

%%%--------------------------
 

figure(1)
legend('N=1','N=2','N=3','N=4' )
grid on
xlabel('Arrival rate of packets (per second)')
ylabel('Network Energy Efficiency (Bit/Joule)')

figure(2)
legend('N=1','N=2','N=3','N=4' )
grid on
xlabel('Arrival rate of packets (per second)')
ylabel('Battery Lifetime (\times reporting period)')

figure(3)
legend('N=1','N=2','N=3','N=4' )
grid on
xlabel('Arrival rate of packets (per second)')
ylabel('Success Rate  ')
%%------------
 

