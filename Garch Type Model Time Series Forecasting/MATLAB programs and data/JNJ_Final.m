%% JNP AR-GARCH-t
logL=0;aict=0;sict=0;
for p=1:5
    for q=1:5
      Mdl = arima(1,0,0); Mdl.Variance=garch(p,q); Mdl.Distribution='T';  
      [EstMdl,EstParamCov,logL,info] = estimate(Mdl,jnj_rtn_is);  
      aict(p,q)=-2*logL+2*(p+q+2); sict(p,q)=-2*logL+log(length(coke_rtn_is))*(p+q+2);
    end
end

% AIC SIF Plot
figure;plot(1:5,aict(1,:),'b+-',1:5,aict(2,:),'k+-',1:5,aict(3,:),'r+-',1:5,aict(4,:),'g+-',1:5,aict(5,:),'y+-');
title('AIC & SIC for Coke AR(1)-GJR-GARCH-t Models');hold on;
plot(1:5,sict(1,:),'b+-',1:5,sict(2,:),'k+-',1:5,sict(3,:),'r+-',1:5,sict(4,:),'g+-',1:5,sict(5,:),'y+-');
legend('AIC,P=1','AIC,P=2','AIC,P=3','AIC,P=4','AIC,P=5','SIC,P=1','SIC,P=2','SIC,P=3','SIC,P=4','SIC,P=5');
xlabel('Number of Q');

% Optimal Model Training
Mdl = arima(1,0,0); Mdl.Variance=garch(1,1); Mdl.Distribution='T'; 
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,jnj_rtn_is);

[jnj_a111t,jnj_v111t,logL] = infer(EstMdl,jnj_rtn_is); jnj_s111t=sqrt(jnj_v111t); %calculate innovations
jnj_e111t=jnj_a111t./jnj_s111t; jnj_e111t_sqr=jnj_e111t.^2;
% G-error transformation
jnj_df=EstMdl.Distribution.DoF;
jnj_e111tg=norminv(tcdf(sqrt(jnj_df)/sqrt(jnj_df-2)*jnj_e111t,jnj_df));

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(jnj_e111tg, [10 15], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(jnj_e111tg.^2, [10 15], 0.05, [5 10])

% JB Test
[skewness(jnj_e111tg) kurtosis(jnj_e111tg)]
[h, p] = jbtest(jnj_e111tg)

%% JNP AR-GJR-GARCH-t
logL=0;aict=0;sict=0;
for p=1:5
    for q=1:5
      Mdl = arima(1,0,0); Mdl.Variance=gjr(p,q); Mdl.Distribution='T';  
      [EstMdl,EstParamCov,logL,info] = estimate(Mdl,jnj_rtn_is);  
      aict(p,q)=-2*LogL+2*(p+q+2); sict(p,q)=-2*LogL+log(length(jnj_rtn_is))*(p+q+2);
    end
end
% AIC SIF Plot
figure;plot(1:5,aict(1,:),'b+-',1:5,aict(2,:),'k+-',1:5,aict(3,:),'r+-',1:5,aict(4,:),'g+-',1:5,aict(5,:),'y+-');
title('AIC & SIC for jnj AR(1)-GJR-GARCH-t Models');hold on;
plot(1:5,sict(1,:),'b+-',1:5,sict(2,:),'k+-',1:5,sict(3,:),'r+-',1:5,sict(4,:),'g+-',1:5,sict(5,:),'y+-');
legend('AIC,P=1','AIC,P=2','AIC,P=3','AIC,P=4','AIC,P=5','SIC,P=1','SIC,P=2','SIC,P=3','SIC,P=4','SIC,P=5');
xlabel('Number of Q');

% Optimal Model Training
Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T'; 
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,jnj_rtn_is);
[jnj_a111g,jnj_v111g,logL] = infer(EstMdl,jnj_rtn_is); jnj_s111g=sqrt(jnj_v111g); %calculate innovations
jnj_e111g=jnj_a111g./jnj_s111g; jnj_e111g_sqr=jnj_e111g.^2;
% G-error transformation
jnj_df=EstMdl.Distribution.DoF;
jnj_e111gg=norminv(tcdf(sqrt(jnj_df)/sqrt(jnj_df-2)*jnj_e111g,jnj_df));
jnj_sigGJt=jnj_s111g; jnj_aGJt=jnj_a111g;

% NIC Part
% NIC CURVE with a
% Get required estimated GJR model coefficients
jnj_a1GJt=cell2mat(EstMdl.Variance.ARCH);jnj_b1GJt=cell2mat(EstMdl.Variance.GARCH);
jnj_g1GJt=cell2mat(EstMdl.Variance.Leverage);
jnj_dfGJt=EstMdl.Variance.Distribution.DoF;  
jnj_a0GJt=EstMdl.Variance.Constant; 

%average volatility in each regime
[jnj_a0GJt/(1-jnj_a1GJt-jnj_g1GJt-jnj_b1GJt)  jnj_a0GJt/(1-jnj_a1GJt-jnj_b1GJt)]

% persistense
[jnj_a1GJt+jnj_g1GJt+jnj_b1GJt jnj_a1GJt+jnj_b1GJt]

% Caculate values for NIC plot
jnj_a=min(jnj_aGJt):0.01:max(jnj_aGJt);   % range of plot
jnj_sg=var(jnj_aGJt);	% sample variance of innovations
jnj_sigt=jnj_a0GJt+jnj_a1GJt*jnj_a.^2+jnj_b1GJt*jnj_sg+(jnj_g1GJt.*(jnj_a<0)).*(jnj_a.^2); % asymmetric curve values
jnj_sigt2=jnj_a0GJt+(jnj_a1GJt+jnj_g1GJt/2)*jnj_a.^2+jnj_b1GJt*jnj_sg;	 % symmetric curve values

% NIC plot for a(t-1)
figure;plot(jnj_a,jnj_sigt);axis([min(jnj_aGJt)-0.5 max(jnj_aGJt)+0.5 0 max(jnj_sigt)+1]);
hold on;plot(jnj_a,jnj_sigt2,'r--');
title('AR(1)-GJR-GARCH(1,1)-t NIC curve for a(t-1)');

% compare precited conditional standard deviation values following positive
% and negative shocks of size 2
jnj_a=-2;jnj_sigtm2=jnj_a0GJt+jnj_a1GJt*jnj_a^2+jnj_b1GJt*jnj_sg+(jnj_g1GJt*(jnj_a<0))*(jnj_a^2); 
jnj_a=2;jnj_sigta2=jnj_a1GJt+jnj_a1GJt*jnj_a^2+jnj_b1GJt*jnj_sg+(jnj_g1GJt*(jnj_a<0))*(jnj_a^2);
[jnj_sigtm2 jnj_sigta2 jnj_sigtm2/jnj_sigta2] 


% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(jnj_e111gg, [12 17], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(jnj_e111gg.^2, [12 17], 0.05, [5 10])

% JB Test
[skewness(jnj_e111gg) kurtosis(jnj_e111gg)]
[h, p] = jbtest(jnj_e111gg)

%% JNJ-IGARCH-N
[PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS]=igarch(jnj_rtn_is,1,1,[],[],0);
SHT=sqrt(HT);
jnj_eHT=jnj_rtn_is./SHT;

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(jnj_eHT, [7 12], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(jnj_eHT.^2, [7 12], 0.05, [5 10])

% JB test
[skewness(jnj_eHT) kurtosis(jnj_eHT)]
[h, p] = jbtest(jnj_eHT)