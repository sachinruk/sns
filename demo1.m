clear all
close all
clc

pi=[0.2 0.8];
N=30; d=100; %20 observations and 100 dimensions
X=randn(N,d);
s=mnrnd(1,pi,d);
w=randn(d,1).*s(:,1);
sigmaq=0.1;
y=X*w+sigmaq*randn(N,1);

[muw, gamma, sigmaw,sigmaq]=vb(y,X,sigmaq^2);