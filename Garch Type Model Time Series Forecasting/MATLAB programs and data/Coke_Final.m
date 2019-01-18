%% 1. Coke AR(1)-GARCH(1,1)-t
logL=0;aict=0;sict=0;
for p=1:5
    for q=1:5
      Mdl = arima(1,0,0); Mdl.Variance=garch(p,q); Mdl.Distribution='T';  
      [EstMdl,EstParamCov,logL,info] = estimate(Mdl,coke_rtn_is);  
      aict(p,q)=-2*logL+2*(p+q+2); sict(p,q)=-2*logL+log(length(coke_rtn_is))*(p+q+2);
    end
end

% AIC SIF Plot
figure;plot(1:5,aict(1,:),'b+-','LineWidth',1.5);hold on;
plot(1:5,aict(2,:),'k+-','LineWidth',1.5);
plot(1:5,aict(3,:),'r+-','LineWidth',1.5);
plot(1:5,aict(4,:),'g+-','LineWidth',1.5);
plot(1:5,aict(5,:),'y+-','LineWidth',1.5);
title('AIC & SIC for Coke AR(1)-GJR-GARCH-t Models');
plot(1:5,sict(1,:),'b+-','LineWidth',1.5);
plot(1:5,sict(2,:),'k+-','LineWidth',1.5);
plot(1:5,sict(3,:),'r+-','LineWidth',1.5);
plot(1:5,sict(4,:),'g+-','LineWidth',1.5);
plot(1:5,sict(5,:),'y+-','LineWidth',1.5);
legend('AIC,P=1','AIC,P=2','AIC,P=3','AIC,P=4','AIC,P=5','SIC,P=1','SIC,P=2','SIC,P=3','SIC,P=4','SIC,P=5');
xlabel('Number of Q');

% Optimal Model Training
Mdl = arima(1,0,0); Mdl.Variance=garch(1,1); Mdl.Distribution='T'; 
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,coke_rtn_is);

[coke_a111t,coke_v111t,logL] = infer(EstMdl,coke_rtn_is); coke_s111t=sqrt(coke_v111t); %calculate innovations
coke_e111t=coke_a111t./coke_s111t; coke_e111t_sqr=coke_e111t.^2;
% G-error transformation
coke_df=EstMdl.Distribution.DoF;
coke_e111tg=norminv(tcdf(sqrt(coke_df)/sqrt(coke_df-2)*coke_e111t,coke_df));
% coke_sigGJt=coke_s111t; coke_aGJt=coke_a111j;
% return in-sample forecast
% coke_rtn_is_111tf=forecast(EstMdl,length(coke_rtn_is),'Y0',coke_rtn_is);
% Plot the innovations, conditional standard deviations and log returns above one another
figure;subplot(3,1,1);plot(date_is(2:end,1),coke_a111t);
% xlim([date_is(2) date_is(end)]);
ylim([min(coke_a111t)-1 max(coke_a111t)+1]);                      % set range of y-axis
title('Coke AR(1)-GJR-GARCH(1,1)-t Innovations');
subplot(3,1,2);plot(date_is(2:end,1),coke_s111t);
% xlim([date_is(2) date_is(end)]);
ylim([0 max(coke_s111t)+1]);                      % set range of y-axis
title('Coke AR(1)-GJR-GARCH(1,1)-t Conditional Standard Deviations');
subplot(3,1,3);plot(date_is(2:end,1),coke_rtn_is);
% xlim([date_is(2) date_is(end)]);
ylim([min(coke_rtn_is)-1 max(coke_rtn_is)+1]);                      % set range of y-axis
title('Coke Log returns');

% Shoulder Plot
figure;plot(date_is(2:end,1),coke_s111t,'b');hold on;plot(date_is(2:end,1),coke_rtn_is,'r');legend('Coke AR(1)-GJR-GARCH(1,1)-t Conditional Std','Coke Log Rtn');
title('Conditional Std for Coke AR(1)-GJR-GARCH(1,1)-t & Log Rtn ');


% ACF ACF Squared plot
% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),coke_e111t,'b'); title('Coke AR(1)-GJR-GARCH(1,1)-t Standardised Residuals');
figure; autocorr(coke_e111t,20); title('Coke AR(1)-GJR-GARCH(1,1)-t Residual ACF Plot');
figure;autocorr(coke_e111t_sqr,20); title('Coke AR(1)-GJR-GARCH(1,1)-t Squared Residual ACF Plot');


% hist qq plot t-dist
figure;subplot(2,1,1);hist(coke_e111t,25);
title('AR(1)-GJR-GARCH(1,1)-t Standardised Residuals');
subplot(2,1,2);qqplot(coke_e111t);
title('QQ plot Coke AR(1)-GJR-GARCH(1,1)-t Standardised Residuals');

% hist qq plot N-dist
figure;subplot(2,1,1);hist(coke_e111tg,25);
title('Coke AR(1)-GJR-GARCH(1,1)-t Normal Inversed Standardised Residuals');
subplot(2,1,2);qqplot(coke_e111tg);
title('Coke AR(1)-GJR-GARCH(1,1)-t Normal Inversed Standardised Residuals');

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(coke_e111tg, [9 14], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(coke_e111tg.^2, [9 14], 0.05, [5 10])

% JB Test
[skewness(coke_e111tg) kurtosis(coke_e111tg)]
[h, p] = jbtest(coke_e111tg)

%% 2. Coke AR(1)-GJR-GARCH(1,1)-t
logL=0;aict=0;sict=0;
for p=1:5
    for q=1:5
      Mdl = arima(1,0,0); Mdl.Variance=gjr(p,q); Mdl.Distribution='T';  
      [EstMdl,EstParamCov,logL,info] = estimate(Mdl,coke_rtn_is);  
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
Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T'; 
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,coke_rtn_is);
[coke_a111g,coke_v111g,logL] = infer(EstMdl,coke_rtn_is); coke_s111g=sqrt(coke_v111g); %calculate innovations
coke_e111g=coke_a111g./coke_s111g; coke_e111g_sqr=coke_e111g.^2;
% G-error transformation
coke_df=EstMdl.Distribution.DoF;
coke_e111gg=norminv(tcdf(sqrt(coke_df)/sqrt(coke_df-2)*coke_e111g,coke_df));
coke_sigGJt=coke_s111g; coke_aGJt=coke_a111g;
% in-sample forecast
% coke_rtn_is_111gf=forecast(EstMdl,1,'Y0',coke_rtn_is);

% Plot the innovations, conditional standard deviations and log returns above one another
figure;subplot(3,1,1);plot(date_is(2:end,1),coke_a111g);
% xlim([date_is(2) date_is(end)]);
ylim([min(coke_a111g)-1 max(coke_a111g)+1]);                      % set range of y-axis
title('Coke AR(1)-GJR-GARCH(1,1)-t Innovations');
subplot(3,1,2);plot(date_is(2:end,1),coke_s111g);
% xlim([date_is(2) date_is(end)]);
ylim([0 max(coke_s111g)+1]);                      % set range of y-axis
title('Coke AR(1)-GJR-GARCH(1,1)-t Conditional Standard Deviations');
subplot(3,1,3);plot(date_is(2:end,1),coke_rtn_is);
% xlim([date_is(2) date_is(end)]);
ylim([min(coke_rtn_is)-1 max(coke_rtn_is)+1]);                      % set range of y-axis
title('Coke Log returns');

% Shoulder Plot
figure;plot(date_is(2:end,1),coke_s111g,'b');hold on;plot(date_is(2:end,1),coke_rtn_is,'r');legend('Coke AR(1)-GJR-GARCH(1,1)-t Conditional Std','Coke Log Rtn');
title('Conditional Standard Deviations for Coke AR(1)-GJR-GARCH(1,1)-t & Log Rtn ');

% ACF ACF Squared plot
% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),coke_e111g,'b'); title('Coke AR(1)-GJR-GARCH(1,1)-t Standardised Residuals');
figure; autocorr(coke_e111g,20); title('Coke AR(1)-GJR-GARCH(1,1)-t Residual ACF Plot');
figure;autocorr(coke_e111g_sqr,20); title('Coke AR(1)-GJR-GARCH(1,1)-t Squared Residual ACF Plot');

% NIC Part
% NIC CURVE with a
% Get required estimated GJR model coefficients
coke_a1GJt=cell2mat(EstMdl.Variance.ARCH);coke_b1GJt=cell2mat(EstMdl.Variance.GARCH);
coke_g1GJt=cell2mat(EstMdl.Variance.Leverage);
coke_dfGJt=EstMdl.Variance.Distribution.DoF;  
coke_a0GJt=EstMdl.Variance.Constant; 

%average volatility in each regime
[coke_a0GJt/(1-coke_a1GJt-coke_g1GJt-coke_b1GJt)  coke_a0GJt/(1-coke_a1GJt-coke_b1GJt)]

% persistense
[coke_a1GJt+coke_g1GJt+coke_b1GJt coke_a1GJt+coke_b1GJt]

% Caculate values for NIC plot
coke_a=min(coke_aGJt):0.01:max(coke_aGJt);   % range of plot
coke_sg=var(coke_aGJt);	% sample variance of innovations
coke_sigt=coke_a0GJt+coke_a1GJt*coke_a.^2+coke_b1GJt*coke_sg+(coke_g1GJt.*(coke_a<0)).*(coke_a.^2); % asymmetric curve values
coke_sigt2=coke_a0GJt+(coke_a1GJt+coke_g1GJt/2)*coke_a.^2+coke_b1GJt*coke_sg;	 % symmetric curve values

% NIC plot for a(t-1)
figure;plot(coke_a,coke_sigt);axis([min(coke_aGJt)-0.5 max(coke_aGJt)+0.5 0 max(coke_sigt)+1]);
hold on;plot(coke_a,coke_sigt2,'r--');
title('AR(1)-GJR-GARCH(1,1)-t NIC curve for a(t-1)');

% compare precited conditional standard deviation values following positive
% and negative shocks of size 2
coke_a=-2;coke_sigtm2=coke_a0GJt+coke_a1GJt*coke_a^2+coke_b1GJt*coke_sg+(coke_g1GJt*(coke_a<0))*(coke_a^2); 
coke_a=2;coke_sigta2=coke_a1GJt+coke_a1GJt*coke_a^2+coke_b1GJt*coke_sg+(coke_g1GJt*(coke_a<0))*(coke_a^2);
[coke_sigtm2 coke_sigta2 coke_sigtm2/coke_sigta2] 

% hist qq plot t-dist
figure;subplot(2,1,1);hist(coke_e111g,25);
title('AR(1)-GJR-GARCH(1,1)-t Standardised Residuals');
subplot(2,1,2);qqplot(coke_e111g);
title('QQ plot Coke AR(1)-GJR-GARCH(1,1)-t Standardised Residuals');

% hist qq plot N-dist
figure;subplot(2,1,1);hist(coke_e111gg,25);
title('Coke AR(1)-GJR-GARCH(1,1)-t Normal Inversed Standardised Residuals');
subplot(2,1,2);qqplot(coke_e111gg);
title('Coke AR(1)-GJR-GARCH(1,1)-t Normal Inversed Standardised Residuals');

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(coke_e111gg, [9 14], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(coke_e111gg.^2, [9 14], 0.05, [5 10])

% JB Test
[skewness(coke_e111gg) kurtosis(coke_e111gg)]
[h, p] = jbtest(coke_e111gg)
%% I-GARCH(1,1)-N
[PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS]=igarch(coke_rtn_is,1,1,[],[],0);
SHT=sqrt(HT);
coke_eHT=coke_rtn_is./SHT;

% Shoulder Plot
figure;plot(date_is(2:end,1),coke_rtn_is,'r');hold on;
plot(date_is(2:end,1),SHT,'b','LineWidth',1.5);
legend('Coke Log Rtn','Coke AR(1)-GJR-GARCH(1,1)-t Conditional Std');
title('Conditional Std for Coke IGARCH-N & Log Rtn ');

% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),coke_eHT,'b'); title('Coke IGARCH Standardised Residuals');
figure; autocorr(coke_eHT,20); title('Coke IGARCH Residual ACF Plot');
figure;autocorr(coke_eHT.^2,20); title('Coke IGARCH Squared Residual ACF Plot');

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(coke_eHT, [7 12], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(coke_eHT.^2, [7 12], 0.05, [5 10])

% hist qq plot t-dist
figure;subplot(2,1,1);hist(coke_eHT,25);
title('Histogram of Coke IGARCH Standardised Residuals');
subplot(2,1,2);qqplot(coke_eHT);
title('QQ plot Coke IGARCH Standardised Residuals');


% JB test
[skewness(coke_eHT) kurtosis(coke_eHT)]
[h, p] = jbtest(coke_eHT)
%% Coke Ad-hoc(25)
csvwrite('coke_rtn_is.csv',coke_rtn_is);
coke_s_hs_is=sqrt(coke_v_hs_is);
%% Model Comparison
% whole c-std estimate
figure;plot(date_is(2:end),coke_rtn_is_111tf,'b');hold on;
plot(date_is(2:end),coke_rtn_is_111gf,'r');
plot(date_is(2:end),SHT,'k');
plot(date_is(26:end),coke_rtn_hs_is,'m');
plot(date_is(2:end),coke_rtn_is,'y');
legend('AR(1)-GARCH(1,1)-t','AR(1)-GJR-GARCH(1,1)-t','I-GARCH(1,1)-N','HS-25','Log Return','Location','southwest');
title('Coke Conditional Std In-Sample Forecast Comparsion'); ylabel('Return');xlabel('Date');hold off;

% the latest 2 year c-std estimate
figure;plot(date_is(end-499:end),coke_rtn_is_111tf(end-499:end),'b');hold on;
plot(date_is(end-499:end),coke_rtn_is_111gf(end-499:end),'r');
plot(date_is(end-499:end),SHT(end-499:end),'k');
plot(date_is(end-499:end),coke_rtn_hs_is(end-499:end),'m');
plot(date_is(end-499:end),coke_rtn_is(end-499:end),'c');
legend('AR(1)-GARCH(1,1)-t','AR(1)-GJR-GARCH(1,1)-t','I-GARCH(1,1)-N','HS-25','Log Return','Location','southwest');
title('Coke Conditional Latest 1y Std In-Sample Forecast Comparsion'); ylabel('Return');xlabel('Date');hold off;

% whole c-return estimate
figure;plot(date_is(2:end),coke_s111t,'b');hold on;
plot(date_is(2:end),coke_s111g,'r');
plot(date_is(2:end),SHT,'k');
plot(date_is(26:end),coke_s_hs_is,'m');
plot(date_is(2:end),coke_rtn_is,'y');
legend('AR(1)-GARCH(1,1)-t','AR(1)-GJR-GARCH(1,1)-t','I-GARCH(1,1)-N','HS-25','Log Return','Location','southwest');
title('Coke Conditional Std In-Sample Forecast Comparsion'); ylabel('Return');xlabel('Date');hold off;

% the latest 2 year c-return estimate
figure;plot(date_is(end-499:end),coke_s111t(end-499:end),'b');hold on;
plot(date_is(end-499:end),coke_s111g(end-499:end),'r');
plot(date_is(end-499:end),SHT(end-499:end),'k');
plot(date_is(end-499:end),coke_s_hs_is(end-499:end),'m');
plot(date_is(end-499:end),coke_rtn_is(end-499:end),'c');
legend('AR(1)-GARCH(1,1)-t','AR(1)-GJR-GARCH(1,1)-t','I-GARCH(1,1)-N','HS-25','Log Return','Location','southwest');
title('Coke Conditional Latest 1y Std In-Sample Forecast Comparsion'); ylabel('Return');xlabel('Date');hold off;