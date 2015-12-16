clear all
close all
clc

rng(1);
pi=[0.2 0.8];
N=50; d=100; %20 observations and 100 dimensions
X=randn(N,d);
s=mnrnd(1,pi,d);
w=randn(d,1).*s(:,1);
sigmaq=0.1;
y=X*w+sigmaq*randn(N,1);

[muw, gamma, Es,sigma2q]=vb3(y,X,sigmaq^2,300);
figure; plot(w);
hold on; plot(muw.*(Es>0.5),'r')
norm(w-X\y)
norm(w-muw.*(Es>0.5))

% [muw, gamma, Es,sigma2q]=vb2(y,X,sigmaq^2, 1000);
% norm(w-muw.*Es)
% figure; plot(w);
% hold on; plot(muw.*Es,'r')
% norm(w-X\y)

% 
% [muw2]=vbRVM(y,X,100);
% norm(w-muw2)
% figure; plot(w);
% hold on; plot(muw2,'m')

% model = slrCreate(X, y,  'Gaussian');
% % run paired variational mean field inference
% iters = 300;  % number of EM iterations
% dispF = 0;    % display lower bound during optimization
% [model vardist, F] = slrPairMF(model, iters, dispF);
% 
% figure; plot(w);
% hold on; plot(vardist.muW,'r')
% norm(w-X\y)
% norm(w-vardist.muW')
% figure; plot(w);
% hold on; plot(vardist.muW.*(vardist.gamma>0.5),'r')
% norm(w-(vardist.muW.*(vardist.gamma>0.5))')
