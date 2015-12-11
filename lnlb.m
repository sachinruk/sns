function F=lnlb(y, X, Es,muw, sigma2qw, sigma2w, sigma2q, alpha)
[N,M]=size(X);
Esw=Es.*muw;
Esw2=Es.*(muw.^2+sigma2qw);
XTX=X'*X;
XTXEsw2=bsxfun(@times,bsxfun(@times,XTX,Esw),Esw');
XEsw=bsxfun(@times,X,Esw');
Ew2=Esw2+(1-Es).*sigma2w;
sumEs=sum(Es);

% error_m=y-2*sum(XEsw,2);
% s_ratio=sigma2q/sigma2w;
% proj_x=sum(X.^2)'+s_ratio;
% sigma2q=(y'*error_m+sum(diag(XTX).*Esw2)+...
%             sum(XTXEsw2(:))-sum(diag(XTXEsw2)))/N;
% 
F0=N*log(2*pi*sigma2q);
F1=(y'*y);
F2=sum(y'*XEsw);
F3=sum(diag(XTX).*Esw2);
F4=sum(XTXEsw2(:))-sum(diag(XTXEsw2));
F5=M*log(2*pi*sigma2w)+sum(Ew2)/sigma2w;
F6=log(alpha)*sumEs+log(1-alpha)*(M-sumEs);

E1=M+M*log(2*pi*sigma2w)-log(sigma2w)*sumEs+sum(Es.*log(sigma2qw));
E2=sum((1-Es).*log(1-Es+(Es==1))+Es.*log(Es+(Es==0)));

F=-0.5*F0+(-0.5*F1+F2-0.5*F3-0.5*F4)/sigma2q-0.5*F5+F6+0.5*E1-E2;
% D=-0.5*(N*log(sigma2q)+(y'*error_m+sum(diag(XTX).*Esw2)+...
%     sum(XTXEsw2(:))-sum(diag(XTXEsw2)))/sigma2q+M*log(sigma2w)+...
%     sum(Ew2)/sigma2w)+log(pi)*sum(Es)+log(1-pi)*sum(1-Es);

% EslnEs=(1-Es).*log(1-Es+(Es==1))+Es.*log(Es+(Es==0));
% E=0.5*(M*log(sigma2w)-log(sigma2w)*sum(Es)+sum(Es.*log(sigma2qw)))...
%     -sum(EslnEs);
% F=D+E;