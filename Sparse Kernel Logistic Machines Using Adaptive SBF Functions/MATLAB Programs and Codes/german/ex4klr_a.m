clear;
gamma=0.005;


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
    eta=0.02; %learning rate;
    eta1=0.5;
    eta0=0.02; % The Learning rate of gradient ascend for model location and scale
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
    a = 0.1*randn(M,1);
    b = 0;
    
    a_tilde = [b; a];
    Xk = [ones(N,1) Phi];  %  pseudo K
    %http://openclassroom.stanford.edu/MainFolder/DocumentPage.php?course=MachineLearning&doc=exercises/ex5/ex5.html
    % https://statacumen.com/teach/SC1/SC1_11_LogisticRegression.pdf
    
    Ix = [0 zeros(1, M); zeros(M,1) eye(M)];
    
    gd_steps = 2;
    for ii = 1:gd_steps
        f = Xk * a_tilde;
        p = 1 ./ (1 + exp(-f));
        V = (1 - p) .* p;
        a_tilde = a_tilde - (Xk' * bsxfun(@times, Xk, V) + gamma * Ix) \ ( Xk' * (p - t) + gamma * [0; a_tilde(2:end)]);
    end
    
    %haty = K * a + b;
    haty = Xk * a_tilde;
    hatp = 1 ./ (1 + exp(-haty));
    %theta2=Phi'*(a.*t);
    tpre= hatp>0.5;
    sum(abs(t-tpre))/ N;   % error rate training set
    
    
    for Iter=1:100
        b = a_tilde(1);
        a = a_tilde(2:end);
        for j=1:M;
            for k=1:N
                id(k)= sum(Mu(j,:).*abs(X(k,:)-C(j,:))) < 1;
                tempC(k,:)=   Mu(j,:).*sign(X(k,:)-C(j,:) )*id(k);   % To be revised according to new phis
                tempMu(k,:)= - abs( C(j,:)- X(k,:) )*id(k);          % To be revised according to new phis
            end;
            for i=1:m
                ParE1(i)=    sign(haty)' * tempC(:,i)*a(j);
                ParE2(i)=    sign(haty)' * tempMu(:,i)*a(j);
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
                Phi(k,j)= max(0, 1-sum(Mu(j,:).*abs(X(k,:)-C(j,:))));  % To be revised according to new phis
            end;
        end
        
        Xk = [ones(N,1) Phi];  %  pseudo K
        for ii = 1:gd_steps
            f = Xk * a_tilde;
            p = 1 ./ (1 + exp(-f));
            V = (1 - p) .* p;
            a_tilde = a_tilde - (Xk' * bsxfun(@times, Xk, V) + gamma * Ix) \ ( Xk' * (p - t) + gamma * [0; a_tilde(2:end)]);
        end
        
        haty = Xk * a_tilde;
        hatp = 1 ./ (1 + exp(-haty));
        %theta2=Phi'*(a.*t);
        tpre= hatp>0.5;
        
        sum(abs(t-tpre))/ N;   % error rate training set
        % J(Iter)=(ones(N,1)-Omega*theta(2:N+1)- theta(2:N+1)/gamma-t*theta(1))'*(ones(N,1)-Omega*theta(2:N+1)- theta(2:N+1)/gamma-t*theta(1));
        % J(Iter)= sum(abs( haty)) ;
        J(Iter)= sum(abs( haty)) ;
    end;
    
    time(ide)=toc
    
    b = a_tilde(1);
    a = a_tilde(2:end);
    
    
    
    
    
    
    
    
    
    
    %  Insert algorithm here
    % test data set;
    N1=length(ttest);
    Xtest=xu;
    
    for i=1:N1
        for j=1:M
            phitest(i,j)=max(0,1-sum(Mu(j,:).*abs(Xtest(i,:)-C(j,:))));
        end;
    end;
    
    hatytest=phitest*a+b;
    prob = exp(hatytest) ./(1 + exp(hatytest));
    tpretest = prob > 0.5;
    Mis(ide)=sum(abs(ttest-tpretest))/N1
end;



