clear;
gamma=0.5;


for ide=1:100;
    la=sprintf('%d',ide);
    traindata=['german_train_data_' la '.asc'];
    xla=load(traindata);
    
    trainlabel=['german_train_labels_' la '.asc'];
    yla=load(trainlabel);
    t=yla; t(yla==-1) = 0;
    
    testdata=['german_test_data_' la '.asc'];
    xu=load(testdata);
    
    testlabel=['german_test_labels_' la '.asc'];
    yu=load(testlabel);
    ttest=yu; ttest(yu==-1) = 0;
    
    ndata=size(xla,1);
    ndataU=size(xu,1);
    
    
    
    
    tic
    
    randn('seed',20);
    rand('seed',20);
    eta=0.002; %learning rate;
    eta1=0.05;
    eta0=0.0002; % The Learning rate of gradient ascend for model location and scale
    N=ndata;  %number of data samples
    M=3;  %number of centers
    m=size(xu,2);   % input dimension
    
    
    mu=0.2;   %initial shaping parameter
    
    
    X= xla;  % Data
    
    
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
    
    % Algorithm 2
    % Initialisation
    alpha   = zeros(N, 1);  % means a
    one     = ones(N, 1);
    K       = Phi * Phi';   %evaluate(net.kernel, x, x);
    eta     = zeros(N, 1);  % Means the model output
    p      = exp(eta)./(1+exp(eta));   %net.invlink(eta);  % Means p
    
    gd_steps = 3;
    for ii = 1:gd_steps
        w     = max(exp(eta)./(1+exp(eta)).^2, 1e-6);
        z     = eta + (t - p)./w;
        R     = chol(K + diag(gamma./w));
        zeta  = R\(R'\one);
        xi    = R\(R'\z);
        delta = 1./sum(zeta);
        b     = sum(xi)*delta;
        alpha = xi - zeta*b;
        eta   = K*alpha + b;
        p    = exp(eta)./(1+exp(eta));
        % code adopted from gkmtool
    end
    
    haty = eta;
    theta2=Phi'*alpha;
    hatp = 1 ./ (1 + exp(-haty));
    tpre= hatp>0.5;
    sum(abs(t-tpre))/ N;   % error rate training set
    
    
    for Iter=1:100
        for j=1:M;
            for k=1:N
                id(k)= sum(Mu(j,:).*abs(X(k,:)-C(j,:))) < 1;
                tempC(k,:)=   Mu(j,:).*sign(X(k,:)-C(j,:) )*id(k);
                tempMu(k,:)= - abs( C(j,:)- X(k,:) )*id(k);
            end;
            for i=1:m
                ParE1(i)=    sign(haty)' * (Phi(:,j)*   tempC(:,i)' +   tempC(:,i)* Phi(:,j)')* (alpha)   ;
                ParE2(i)=    sign(haty)' * (Phi(:,j)*   tempMu(:,i)'+   tempMu(:,i)*Phi(:,j)')* (alpha);
            end;
            if norm( ParE1) ~=0
                C(j,:)= C(j,:)+eta0*ParE1 / norm( ParE1) ;
            end;
            if norm( ParE2) ~=0
                Mu(j,:)=Mu(j,:)+eta0*ParE2 / norm( ParE2) ;
            end;
            for it=1:m;
                Mu(j,it)  =max(0, Mu(j,it));
            end;
            for k=1:N
                Phi(k,j)= max(0, 1-sum(Mu(j,:).*abs(X(k,:)-C(j,:))));
            end;
        end
        
        K       = Phi * Phi';   %evaluate(net.kernel, x, x);
        p      = exp(eta)./(1+exp(eta));   %net.invlink(eta);  % Means p
        for ii = 1:gd_steps
            w     = max(exp(eta)./(1+exp(eta)).^2, 1e-6);
            z     = eta + (t - p)./w;
            R     = chol(K + diag(gamma./w));
            zeta  = R\(R'\one);
            xi    = R\(R'\z);
            delta = 1./sum(zeta);
            b     = sum(xi)*delta;
            alpha = xi - zeta*b;
            eta   = K*alpha + b;
            p    = exp(eta)./(1+exp(eta));
            % code adopted from gkmtool
        end
        theta2=Phi'*alpha;
        haty = eta;
        hatp = 1 ./ (1 + exp(-haty));
        tpre= hatp>0.5;
        
        sum(abs(t-tpre))/ N;   % error rate training set
        % J(Iter)=(ones(N,1)-Omega*theta(2:N+1)- theta(2:N+1)/gamma-t*theta(1))'*(ones(N,1)-Omega*theta(2:N+1)- theta(2:N+1)/gamma-t*theta(1));
        % J(Iter)= sum(abs( haty)) ;
        J(Iter)= sum(abs( haty)) ;
    end;
    time(ide)=toc
    
    
    %  Insert algorithm here
    % test data set;
    N1=length(ttest);
    Xtest=xu;
    
    for i=1:N1
        for j=1:M
            phitest(i,j)=max(0,1-sum(Mu(j,:).*abs(Xtest(i,:)-C(j,:))));
        end;
    end;
    
    hatytest=phitest*theta2+b;
    prob = exp(hatytest) ./(1 + exp(hatytest));
    tpretest = prob > 0.5;
    Mis(ide)=sum(abs(ttest-tpretest))/N1   % error rate test set
end;



