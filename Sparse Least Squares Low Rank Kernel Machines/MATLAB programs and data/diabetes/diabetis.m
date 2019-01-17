clear;
gamma=50000;


for ide=1:100;
    la=sprintf('%d',ide);
    traindata=['diabetis_train_data_' la '.asc'];
    xla=load(traindata);
    
    trainlabel=['diabetis_train_labels_' la '.asc'];
    yla=load(trainlabel);
    t=yla;
    
    testdata=['diabetis_test_data_' la '.asc'];
    xu=load(testdata);
    
    testlabel=['diabetis_test_labels_' la '.asc'];
    yu=load(testlabel);
    ttest=yu;
    
    ndata=size(xla,1);
    ndataU=size(xu,1);
    
    
    
    
    
    randn('seed',20);
    rand('seed',20);
    eta=0.005; %learning rate;
    eta1=0.05;
    N=ndata;  %number of data samples
    M=5;  %number of centers
    m=size(xu,2);   % input dimension
    
    
    mu=0.05;   %initial shaping parameter
    
    
    X= xla;  % Data
    t= yla; %Data label
    
    
    Cind=randperm(N);  %
    C=X(Cind(1:M),:);   % randomly choose initial centers from data.
    
    for iter=1:20000;   %repeatly sampling data
        temp=round(rand(1)*(N-1))+1;
        
        xdraw = X(temp,:) ; %randomly draw data sample;
        
        %Find the nearest center;
        for i=1:M;
            Dis(i)=norm(xdraw-C(i,:));
        end;
        [Cmin,Id]=min(Dis);
        %update C(Id,:)
        C(Id,:)= C(Id,:)+ eta1*(xdraw-C(Id,:));
        
    end;
    
    % Change to kmeans
    
    Mu= mu*ones(M,m); %initial shaping parameter
    
    
    
    
    
    
    
    for k=1:N
        for j=1:M
            Phi(k,j)= max(0,1-sum(Mu(j,:).*abs(X(k,:)-C(j,:))));
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
    sum(abs(t-tpre)/2)/N;  % error rate training set
    
    
    
    for Iter=1:500
        for j=1:M;
            
            for k=1:N
                id(k)= sum(Mu(j,:).*abs(X(k,:)-C(j,:))) < 1;
                tempC(k,:)=   Mu(j,:).*sign(X(k,:)-C(j,:) )*id(k);
                tempMu(k,:)= - abs( C(j,:)- X(k,:) )*id(k);
            end;
            for i=1:m
                ParE1(i)=    (sign(haty)'*    Phi(:,j))*  ( tempC(:,i)'*(t.* theta(2:N+1)))...
                    +  (sign(haty)'*  tempC(:,i))* (Phi(:,j)' * (t.* theta(2:N+1)));
                
                ParE2(i)=    ( sign(haty)'  *  Phi(:,j))*   (tempMu(:,i)'*(t.* theta(2:N+1)))...
                    +  (sign(haty)'* tempMu(:,i))*(Phi(:,j)' * (t.* theta(2:N+1)));
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
                Phi(k,j)= max(0, 1-sum(Mu(j,:).*abs(X(k,:)-C(j,:))));
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
        
        theta=q- P*(tildePhi*((inv(eye(M)+tildePhi'*P*tildePhi)*(tildePhi'*q))));
        
        theta2=Phi'*(theta(2:N+1).*t);
        
        haty =Phi*theta2+theta(1);
        tpre=sign(haty);
        sum(abs(t-tpre)/2)/N;   % error rate training set
        
        % J(Iter)=(ones(N,1)-Omega*theta(2:N+1)- theta(2:N+1)/gamma-t*theta(1))'*(ones(N,1)-Omega*theta(2:N+1)- theta(2:N+1)/gamma-t*theta(1));
        % J(Iter)= sum(abs( haty)) ;
        J(Iter)= sum(abs( haty)) ;
    end;
    
    
    
    
    
    
    
    
    
    
    %  Insert algorithm here
    % test data set;
    
    for k=1:ndataU
        for j=1:M
            PhiU(k,j)= max(0,1-sum(Mu(j,:).*abs(xu(k,:)-C(j,:))));
        end;
    end;
    
    
    yU=PhiU*theta2+theta(1);
    
    tildeUT= sign(yU);
    Mis(ide)=sum(abs(ttest-tildeUT)/2)/ndataU
    
end;



