%% 1. MSFT AR(1)-GARCH(1,1)-t
logL=0;aict=0;sict=0;
for p=1:5
    for q=1:5
      Mdl = arima(1,0,0); Mdl.Variance=garch(p,q); Mdl.Distribution='T';  
      [EstMdl,EstParamCov,logL,info] = estimate(Mdl,msft_rtn_is);  
      aict(p,q)=-2*logL+2*(p+q+2); sict(p,q)=-2*logL+log(length(msft_rtn_is))*(p+q+2);
    end
end

% AIC SIF Plot
figure;plot(1:5,aict(1,:),'b+-',1:5,aict(2,:),'k+-',1:5,aict(3,:),'r+-',1:5,aict(4,:),'g+-',1:5,aict(5,:),'y+-');
title('AIC & SIC for MSFT AR(1)-GJR-GARCH-t Models');hold on;
plot(1:5,sict(1,:),'b+-',1:5,sict(2,:),'k+-',1:5,sict(3,:),'r+-',1:5,sict(4,:),'g+-',1:5,sict(5,:),'y+-');
legend('AIC,P=1','AIC,P=2','AIC,P=3','AIC,P=4','AIC,P=5','SIC,P=1','SIC,P=2','SIC,P=3','SIC,P=4','SIC,P=5');
xlabel('Number of Q');

% Optimal Model Training
Mdl = arima(1,0,0); Mdl.Variance=garch(1,1); Mdl.Distribution='T'; 
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,msft_rtn_is);

[msft_a111t,msft_v111t,logL] = infer(EstMdl,msft_rtn_is); msft_s111t=sqrt(msft_v111t); %calculate innovations
msft_e111t=msft_a111t./msft_s111t; msft_e111t_sqr=msft_e111t.^2;
% G-error transformation
msft_df=EstMdl.Distribution.DoF;
msft_e111tg=norminv(tcdf(sqrt(msft_df)/sqrt(msft_df-2)*msft_e111t,msft_df));
% msft_sigGJt=msft_s111t; msft_aGJt=msft_a111j;
% return in-sample forecast
% msft_rtn_is_111tf=forecast(EstMdl,length(msft_rtn_is),'Y0',msft_rtn_is);
% Plot the innovations, conditional standard deviations and log returns above one another
figure;subplot(3,1,1);plot(date_is(2:end,1),msft_a111t);
% xlim([date_is(2) date_is(end)]);
ylim([min(msft_a111t)-1 max(msft_a111t)+1]);                      % set range of y-axis
title('msft AR(1)-GJR-GARCH(1,1)-t Innovations');
subplot(3,1,2);plot(date_is(2:end,1),msft_s111t);
% xlim([date_is(2) date_is(end)]);
ylim([0 max(msft_s111t)+1]);                      % set range of y-axis
title('msft AR(1)-GJR-GARCH(1,1)-t Conditional Standard Deviations');
subplot(3,1,3);plot(date_is(2:end,1),msft_rtn_is);
% xlim([date_is(2) date_is(end)]);
ylim([min(msft_rtn_is)-1 max(msft_rtn_is)+1]);                      % set range of y-axis
title('msft Log returns');

% Shoulder Plot
figure;plot(date_is(2:end,1),msft_s111t,'b');hold on;plot(date_is(2:end,1),msft_rtn_is,'r');legend('msft AR(1)-GJR-GARCH(1,1)-t Conditional Std','msft Log Rtn');
title('Conditional Std for msft AR(1)-GJR-GARCH(1,1)-t & Log Rtn ');


% ACF ACF Squared plot
% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),msft_e111t,'b'); title('msft AR(1)-GJR-GARCH(1,1)-t Standardised Residuals');
figure; autocorr(msft_e111t,20); title('msft AR(1)-GJR-GARCH(1,1)-t Residual ACF Plot');
figure;autocorr(msft_e111t_sqr,20); title('msft AR(1)-GJR-GARCH(1,1)-t Squared Residual ACF Plot');


% hist qq plot t-dist
figure;subplot(2,1,1);hist(msft_e111t,25);
title('AR(1)-GJR-GARCH(1,1)-t Standardised Residuals');
subplot(2,1,2);qqplot(msft_e111t);
title('QQ plot msft AR(1)-GJR-GARCH(1,1)-t Standardised Residuals');

% hist qq plot N-dist
figure;subplot(2,1,1);hist(msft_e111tg,25);
title('msft AR(1)-GJR-GARCH(1,1)-t Normal Inversed Standardised Residuals');
subplot(2,1,2);qqplot(msft_e111tg);
title('msft AR(1)-GJR-GARCH(1,1)-t Normal Inversed Standardised Residuals');

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(msft_e111tg, [9 14], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(msft_e111tg.^2, [9 14], 0.05, [5 10])

% JB Test
[skewness(msft_e111tg) kurtosis(msft_e111tg)]
[h, p] = jbtest(msft_e111tg)

%% 2. MSFT AR(1)-GJR-GARCH(1,1)-t
logL=0;aict=0;sict=0;
for p=1:5
    for q=1:5
      Mdl = arima(1,0,0); Mdl.Variance=gjr(p,q); Mdl.Distribution='T';  
      [EstMdl,EstParamCov,logL,info] = estimate(Mdl,msft_rtn_is);  
      aict(p,q)=-2*logL+2*(p+q+2); sict(p,q)=-2*logL+log(length(msft_rtn_is))*(p+q+2);
    end
end
% AIC SIF Plot
figure;plot(1:5,aict(1,:),'b+-',1:5,aict(2,:),'k+-',1:5,aict(3,:),'r+-',1:5,aict(4,:),'g+-',1:5,aict(5,:),'y+-');
title('AIC & SIC for MSFT AR(1)-GJR-GARCH-t Models');hold on;
plot(1:5,sict(1,:),'b+-',1:5,sict(2,:),'k+-',1:5,sict(3,:),'r+-',1:5,sict(4,:),'g+-',1:5,sict(5,:),'y+-');
legend('AIC,P=1','AIC,P=2','AIC,P=3','AIC,P=4','AIC,P=5','SIC,P=1','SIC,P=2','SIC,P=3','SIC,P=4','SIC,P=5');
xlabel('Number of Q');

% Optimal Model Training
Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T'; 
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,msft_rtn_is);
[msft_a111g,msft_v111g,logL] = infer(EstMdl,msft_rtn_is); msft_s111g=sqrt(msft_v111g); %calculate innovations
msft_e111g=msft_a111g./msft_s111g; msft_e111g_sqr=msft_e111g.^2;
% G-error transformation
msft_df=EstMdl.Distribution.DoF;
msft_e111gg=norminv(tcdf(sqrt(msft_df)/sqrt(msft_df-2)*msft_e111g,msft_df));
msft_sigGJt=msft_s111g; msft_aGJt=msft_a111g;
% in-sample forecast
% msft_rtn_is_111gf=forecast(EstMdl,1,'Y0',msft_rtn_is);

% Plot the innovations, conditional standard deviations and log returns above one another
figure;subplot(3,1,1);plot(date_is(2:end,1),msft_a111g);
% xlim([date_is(2) date_is(end)]);
ylim([min(msft_a111g)-1 max(msft_a111g)+1]);                      % set range of y-axis
title('msft AR(1)-GJR-GARCH(1,1)-t Innovations');
subplot(3,1,2);plot(date_is(2:end,1),msft_s111g);
% xlim([date_is(2) date_is(end)]);
ylim([0 max(msft_s111g)+1]);                      % set range of y-axis
title('msft AR(1)-GJR-GARCH(1,1)-t Conditional Standard Deviations');
subplot(3,1,3);plot(date_is(2:end,1),msft_rtn_is);
% xlim([date_is(2) date_is(end)]);
ylim([min(msft_rtn_is)-1 max(msft_rtn_is)+1]);                      % set range of y-axis
title('msft Log returns');

% Shoulder Plot
figure;plot(date_is(2:end,1),msft_s111g,'b');hold on;plot(date_is(2:end,1),msft_rtn_is,'r');legend('msft AR(1)-GJR-GARCH(1,1)-t Conditional Std','msft Log Rtn');
title('Conditional Standard Deviations for msft AR(1)-GJR-GARCH(1,1)-t & Log Rtn ');

% ACF ACF Squared plot
% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),msft_e111g,'b'); title('msft AR(1)-GJR-GARCH(1,1)-t Standardised Residuals');
figure; autocorr(msft_e111g,20); title('msft AR(1)-GJR-GARCH(1,1)-t Residual ACF Plot');
figure;autocorr(msft_e111g_sqr,20); title('msft AR(1)-GJR-GARCH(1,1)-t Squared Residual ACF Plot');

% NIC Part
% NIC CURVE with a
% Get required estimated GJR model coefficients
msft_a1GJt=cell2mat(EstMdl.Variance.ARCH);msft_b1GJt=cell2mat(EstMdl.Variance.GARCH);
msft_g1GJt=cell2mat(EstMdl.Variance.Leverage);
msft_dfGJt=EstMdl.Variance.Distribution.DoF;  
msft_a0GJt=EstMdl.Variance.Constant; 

%average volatility in each regime
[msft_a0GJt/(1-msft_a1GJt-msft_g1GJt-msft_b1GJt)  msft_a0GJt/(1-msft_a1GJt-msft_b1GJt)]

% persistense
[msft_a1GJt+msft_g1GJt+msft_b1GJt msft_a1GJt+msft_b1GJt]

% Caculate values for NIC plot
msft_a=min(msft_aGJt):0.01:max(msft_aGJt);   % range of plot
msft_sg=var(msft_aGJt);	% sample variance of innovations
msft_sigt=msft_a0GJt+msft_a1GJt*msft_a.^2+msft_b1GJt*msft_sg+(msft_g1GJt.*(msft_a<0)).*(msft_a.^2); % asymmetric curve values
msft_sigt2=msft_a0GJt+(msft_a1GJt+msft_g1GJt/2)*msft_a.^2+msft_b1GJt*msft_sg;	 % symmetric curve values

% NIC plot for a(t-1)
figure;plot(msft_a,msft_sigt);axis([min(msft_aGJt)-0.5 max(msft_aGJt)+0.5 0 max(msft_sigt)+1]);
hold on;plot(msft_a,msft_sigt2,'r--');
title('AR(1)-GJR-GARCH(1,1)-t NIC curve for a(t-1)');

% compare precited conditional standard deviation values following positive
% and negative shocks of size 2
msft_a=-2;msft_sigtm2=msft_a0GJt+msft_a1GJt*msft_a^2+msft_b1GJt*msft_sg+(msft_g1GJt*(msft_a<0))*(msft_a^2); 
msft_a=2;msft_sigta2=msft_a1GJt+msft_a1GJt*msft_a^2+msft_b1GJt*msft_sg+(msft_g1GJt*(msft_a<0))*(msft_a^2);
[msft_sigtm2 msft_sigta2 msft_sigtm2/msft_sigta2] 

% hist qq plot t-dist
figure;subplot(2,1,1);hist(msft_e111g,25);
title('AR(1)-GJR-GARCH(1,1)-t Standardised Residuals');
subplot(2,1,2);qqplot(msft_e111g);
title('QQ plot msft AR(1)-GJR-GARCH(1,1)-t Standardised Residuals');

% hist qq plot N-dist
figure;subplot(2,1,1);hist(msft_e111gg,25);
title('msft AR(1)-GJR-GARCH(1,1)-t Normal Inversed Standardised Residuals');
subplot(2,1,2);qqplot(msft_e111gg);
title('msft AR(1)-GJR-GARCH(1,1)-t Normal Inversed Standardised Residuals');

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(msft_e111gg, [9 14], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(msft_e111gg.^2, [9 14], 0.05, [5 10])

% JB Test
[skewness(msft_e111gg) kurtosis(msft_e111gg)]
[h, p] = jbtest(msft_e111gg)
%% MSFT-IGARCH-N
[PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS]=igarch(msft_rtn_is,1,1,[],[],0);
SHT=sqrt(HT);
msft_eHT=msft_rtn_is./SHT;

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(msft_eHT, [7 12], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(msft_eHT.^2, [7 12], 0.05, [5 10])

% JB test
[skewness(msft_eHT) kurtosis(msft_eHT)]
[h, p] = jbtest(msft_eHT)
%% msft Ad-hoc(25)
csvwrite('msft_rtn_is.csv',msft_rtn_is);
msft_s_hs_is=sqrt(msft_v_hs_is);
%% Model Comparison
% whole c-std estimate
figure;plot(date_is(2:end),msft_rtn_is_111tf,'b');hold on;
plot(date_is(2:end),msft_rtn_is_111gf,'r');
plot(date_is(2:end),SHT,'k');
plot(date_is(26:end),msft_rtn_hs_is,'m');
plot(date_is(2:end),msft_rtn_is,'y');
legend('AR(1)-GARCH(1,1)-t','AR(1)-GJR-GARCH(1,1)-t','I-GARCH(1,1)-N','HS-25','Log Return','Location','southwest');
title('msft Conditional Std In-Sample Forecast Comparsion'); ylabel('Return');xlabel('Date');hold off;

% the latest 2 year c-std estimate
figure;plot(date_is(end-499:end),msft_rtn_is_111tf(end-499:end),'b');hold on;
plot(date_is(end-499:end),msft_rtn_is_111gf(end-499:end),'r');
plot(date_is(end-499:end),SHT(end-499:end),'k');
plot(date_is(end-499:end),msft_rtn_hs_is(end-499:end),'m');
plot(date_is(end-499:end),msft_rtn_is(end-499:end),'c');
legend('AR(1)-GARCH(1,1)-t','AR(1)-GJR-GARCH(1,1)-t','I-GARCH(1,1)-N','HS-25','Log Return','Location','southwest');
title('msft Conditional Latest 1y Std In-Sample Forecast Comparsion'); ylabel('Return');xlabel('Date');hold off;

% whole c-return estimate
figure;plot(date_is(2:end),msft_s111t,'b');hold on;
plot(date_is(2:end),msft_s111g,'r');
plot(date_is(2:end),SHT,'k');
plot(date_is(26:end),msft_s_hs_is,'m');
plot(date_is(2:end),msft_rtn_is,'y');
legend('AR(1)-GARCH(1,1)-t','AR(1)-GJR-GARCH(1,1)-t','I-GARCH(1,1)-N','HS-25','Log Return','Location','southwest');
title('msft Conditional Std In-Sample Forecast Comparsion'); ylabel('Return');xlabel('Date');hold off;

% the latest 2 year c-return estimate
figure;plot(date_is(end-499:end),msft_s111t(end-499:end),'b');hold on;
plot(date_is(end-499:end),msft_s111g(end-499:end),'r');
plot(date_is(end-499:end),SHT(end-499:end),'k');
plot(date_is(end-499:end),msft_s_hs_is(end-499:end),'m');
plot(date_is(end-499:end),msft_rtn_is(end-499:end),'c');
legend('AR(1)-GARCH(1,1)-t','AR(1)-GJR-GARCH(1,1)-t','I-GARCH(1,1)-N','HS-25','Log Return','Location','southwest');
title('msft Conditional Latest 1y Std In-Sample Forecast Comparsion'); ylabel('Return');xlabel('Date');hold off;