clear;
load sytr.dat
load synte.dat;
temp1=find(sytr(:,3)==1);
temp2=find(sytr(:,3)==-1);
figure(2);hold on;
plot(sytr(temp2,1),sytr(temp2,2),'ro');
plot(sytr(temp1,1),sytr(temp1,2),'kx');
X=sytr(:,1:2);
t=sytr(:,3);
N=length(t);
randn('seed',130);
rand('seed',23);
M=3;

eta=0.0008;
% eta1=0.02;
m=2;
mu=0.2;  %initial shaping parameter

gamma=150;

% Cind=randperm(N);  %
% C=X(Cind(1:M),:);   % randomly choose initial centers from data.
% 
% for iter=1:20000;   %repeatly sampling data
%     temp=round(rand(1)*(N-1))+1;
%     
%     xdraw = X(temp,:) ; %randomly draw data sample;
%     
%     %Find the nearest center;
%     for i=1:M;
%         Dis(i)=norm(xdraw-C(i,:));
%     end;
%     [Cmin,Id]=min(Dis);
%     %update C(Id,:)
%     C(Id,:)= C(Id,:)+ eta1*(xdraw-C(Id,:));
%     
% end;
[idx,C] = kmedoids(X,M); % K-Median Clustering
% [idx,C] = kmeans(X,M); % K-Means Clustering


Mu= mu*ones(M,m); %initial shaping parameter




for k=1:N
    for j=1:M
        Phi(k,j)= exp(-sum(Mu(j,:).*abs(X(k,:)-C(j,:))));
    end;
end;

%Algorithm 2
P=[-1/gamma   t'; t gamma*(N*eye(N)-t*t')]/N;
q=[sum(t); gamma*(N*ones(N,1)-sum(t)*t)]/N;

tildePhi= [];
for j=1:M;
    tildePhi=[tildePhi t.*Phi(:,j)];
end;
Omega=tildePhi*tildePhi';
tildePhi=[zeros(1,M);tildePhi];

theta=q- P*tildePhi*inv(eye(M)+tildePhi'*P*tildePhi)*tildePhi'*q;

theta2=Phi'*(theta(2:N+1).*t);

haty =Phi*theta2+theta(1);
tpre=sign(haty);
sum(abs(t-tpre)/2)/N   % error rate training set



for Iter=1:100
    for j=1:M;
        
        for k=1:N
            tempC(k,:)=   Mu(j,:).*sign(X(k,:)-C(j,:)).*exp(- abs( C(j,:)- X(k,:)).*Mu(j,:));
            tempMu(k,:)= - abs( C(j,:)- X(k,:)).*exp(- abs(C(j,:)- X(k,:)).*Mu(j,:));
        end;
        
        
        for i=1:m
            ParE1(i)=    t'*    (Phi(:,j)*   tempC(:,i)' +   tempC(:,i)* Phi(:,j)')* (t.* theta(2:N+1))   ;
            ParE2(i)=     t'  *( Phi(:,j)*   tempMu(:,i)'+   tempMu(:,i)*Phi(:,j)')* (t.* theta(2:N+1));
        end;
        
        
        if norm( ParE1) ~=0
            C(j,:)= C(j,:)+eta*ParE1 / norm( ParE1) ;
        end;
        if norm( ParE2) ~=0
            Mu(j,:)=Mu(j,:)+eta*ParE2 / norm( ParE2) ;
        end;
        for it=1:m;
            Mu(j,it)  =max(0, Mu(j,it));
        end;
        for k=1:N
            Phi(k,j)= exp(-sum(Mu(j,:).*abs(X(k,:)-C(j,:))));
        end;
        
        
        %  J(Iter)=   ( t.* theta(2:N+1))' *Phi*Phi'*   ( t.* theta(2:N+1))  ;
        
        
        %Algorithm 2
        % P=[-1/gamma   t'; t gamma*(N*eye(N)-t*t')]/N;
        % q=[sum(t); gamma*(N*ones(N,1)-sum(t)*t)]/N;
        
        % tildePhi= [];
        % for j=1:M;
        %   tildePhi=[tildePhi t.*Phi(:,j)];
        % end;
        %
        % tildePhi=[zeros(1,M);tildePhi];
        
    end
    tildePhi= [];
    for j=1:M;
        tildePhi=[tildePhi t.*Phi(:,j)];
    end;
    Omega=tildePhi*tildePhi';
    tildePhi=[zeros(1,M);tildePhi];
    
    theta=q- P*tildePhi*inv(eye(M)+tildePhi'*P*tildePhi)*tildePhi'*q;
    
    theta2=Phi'*(theta(2:N+1).*t);
    
    haty =Phi*theta2+theta(1);
    tpre=sign(haty);
    sum(abs(t-tpre)/2)/N   % error rate training set
    
    % J(Iter)=(ones(N,1)-Omega*theta(2:N+1)- theta(2:N+1)/gamma-t*theta(1))'*(ones(N,1)-Omega*theta(2:N+1)- theta(2:N+1)/gamma-t*theta(1));
    % J(Iter)= sum(abs( haty)) ;
    J(Iter)= sum(abs( haty)) ;
end;
% % generate the decision boundary.

for i1=1:100;
    for j1=1:120;
        xn(i1,j1)=i1*1/50-1 ;
        yn(i1,j1)=j1*1/100;
        tem=[xn(i1,j1) yn(i1,j1)];
        for j=1:M
            phi(j)= exp(-sum(Mu(j,:).*abs(tem-C(j,:)))); 
        end;
        zn(i1,j1)=phi*theta2+theta(1);
    end;
end;
figure(2);
contour(xn,yn,zn, [0 0],'k.','linewidth',2);
xlabel('x_1');ylabel('x_2');
xlim([-1.5 1]); ylim([-0.2 1.2]);
% print('out', '-dpdf');


%test data
ttest=synte(:,3);
N1=length(ttest);
Xtest=synte(:,1:2);
for i=1:N1
    for j=1:M
        
        phitest(i,j)= exp(-sum(Mu(j,:).*abs(Xtest(i,:)-C(j,:)))); 
        
    end;
end;

hatytest=phitest*theta2+theta(1);
tpretest=sign(hatytest);
sum(abs(ttest-tpretest)/2)/N1   % error rate test set


