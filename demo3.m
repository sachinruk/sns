clear all
close all
clc

rng(1);
pi=[0.2 0.8];
N=300; d=100; %20 observations and 100 dimensions
X=randn(N,d);
s=mnrnd(1,pi,d);
w=randn(d,1).*s(:,1);
sigmaq=0.1;
y=X*w+sigmaq*randn(N,1);
z=ones(size(y));
z(y<0)=0;

iter=300;
[muw, sigma2qw, Es]=classify(z,X, iter);
muw2=RVMclassification(z,X, iter);
muw3=ProbitClassification(z,X, iter);
% [muw2,Sigma,A]=classify_RVM(z,X,iter);
norm(w-muw)
norm(w-muw2)
norm(w-muw3)
% norm(z-Es)
% norm(z-A)