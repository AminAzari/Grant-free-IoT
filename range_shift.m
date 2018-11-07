%IHN
clc
clear all
close all

L=45;
Tb=0.001;
Fs=8000;
t=[1/Fs:1/Fs:Tb];
i=0;
fd_r=[-1000:5:1000];

for fd=fd_r
    i=i+1;
sh(i)=FunDrN(L,Tb,Fs,fd);
end

plot(fd_r,sh/L)