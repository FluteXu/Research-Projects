save assignment.mat
% 
% coke_rtn=100*diff(log(COKE(:,6)));
% jnj_rtn=100*diff(log(JNJ(:,6)));
% mmm_rtn=100*diff(log(MMM(:,6)));
% msft_rtn=100*diff(log(MSFT(:,6)));
% ptr_rtn=100*diff(log(PTR(:,6)));
%% data processing
% train test split
% Coke
% coke_rtn_is=coke_rtn(1:end-840);
% coke_rtn_os=coke_rtn(end-839:end);
% % JNJ
% jnj_rtn_is=jnj_rtn(1:end-840);
% jnj_rtn_os=jnj_rtn(end-839:end);
% % MMM
% mmm_rtn_is=mmm_rtn(1:end-840);
% mmm_rtn_os=mmm_rtn(end-839:end);
% % MSFT
% msft_rtn_is=msft_rtn(1:end-840);
% msft_rtn_os=msft_rtn(end-839:end);
% % PTR
% ptr_rtn_is=ptr_rtn(1:end-840);
% ptr_rtn_os=ptr_rtn(end-839:end);
% 
% date_is=dates(1:end-840);
% date_os=dates(end-839:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PTR Model Fitting %%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. AR(1)-GARCH-t
% Firstly consider AR(1)-GARCH(p,q) models with gaussian and student t errors
logL=0;aict=0;sict=0;
for p=1:5
  for q=1:5
% NOTE that in MATLAB p is number of lagged conditional std dev's, q is number of lagged squared innovations
% this notation is different to lectures (p and q are switched) and most of the GARCH literature
%     Mdl = arima('ARLags',1,'Variance',garch(p,q));Mdl.Distribution='Gaussian';
%     [EstMdl,EstParamCov,logL,info] = estimate(Mdl,BHPr,'display','off'); 
%     aic(p,q)=-2*logL+2*(p+q+1);sic1(p,q)=-2*logL+log(length(BHPr))*(p+q+1);
    Mdl = arima('ARLags',1,'Variance',garch(p,q),'Distribution','T');
    [EstMdl,EstParamCov,logL,info] = estimate(Mdl,ptr_rtn_is,'display','off');
    aict(p,q)=-2*logL+2*(p+q+2); sict(p,q)=-2*logL+log(length(ptr_rtn_is))*(p+q+2);
  end
end 
% AIC select AR(1)-GARCH(2,2)-t SIC select AR(1)-GARCH(2,2)-t
%
% AIC&SIC plot
figure;plot(aict,'b+-');title('AIC & SIC for PTR GARCH-t models');hold on;plot(sict,'r+-');
legend('AIC','SIC');

[aict_min,aict_min_index]=min(aict);
[sict_min,sict_min_index]=min(sict);

% optimal AR-GARCH-t model fit
Mdl = arima('ARLags',1,'Variance',garch(2,2),'Distribution','T'); %specify model
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,ptr_rtn_is);
[ptr_a22t,ptr_v22t,logL] = infer(EstMdl,ptr_rtn_is);ptr_s22t=sqrt(ptr_v22t); %calculate innovations
ptr_e22t=ptr_a22t./ptr_s22t; ptr_e22t_sqr=ptr_e22t.^2;
% G-error transformation
ptr_df=EstMdl.Distribution.DoF;
ptr_e22tg=norminv(tcdf(sqrt(ptr_df)/sqrt(ptr_df-2)*ptr_e22t,ptr_df));

% Plot the innovations, conditional standard deviations and log returns above one another
figure;subplot(3,1,1);plot(date_is(2:end,1),ptr_a22t);
% xlim([date_is(2) date_is(end)]);
ylim([min(ptr_a22t)-1 max(ptr_a22t)+1]);                      % set range of y-axis
title('PTR AR-GARCH(2,2)-t Innovations');
subplot(3,1,2);plot(date_is(2:end,1),ptr_s22t);
% xlim([date_is(2) date_is(end)]);
ylim([0 max(ptr_s22t)+1]);                      % set range of y-axis
title('PTR AR-GARCH(2,2)-t Conditional Standard Deviations');
subplot(3,1,3);plot(date_is(2:end,1),ptr_rtn_is);
% xlim([date_is(2) date_is(end)]);
ylim([min(ptr_rtn_is)-1 max(ptr_rtn_is)+1]);                      % set range of y-axis
title('PTR Log returns');

% AR-GARCH(2,2)-T Residual Diagnostics
% Cstd & log rtn
figure;plot(date_is(2:end,1),ptr_s22t,'b');hold on;plot(date_is(2:end,1),ptr_rtn_is,'r');legend('PTR AR-GARCH(2,2)-t Conditional Std','PTR Log Rtn');
title('Conditional Standard Deviations for PTR AR-GARCH(2,2)-t & Log Rtn ');

% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),ptr_e22t,'b'); title('PTR Standardised Residuals');
figure; autocorr(ptr_e22t,20); title('PTR Residual ACF Plot');
figure;autocorr(ptr_e22t_sqr,20); title('PTR Squared Residual ACF Plot');

% hist qq plot t-dist
figure;subplot(2,1,1);hist(ptr_e22t,25);
title('Histogram of PTR AR-GARCH(2,2)-t Standardised Residuals');
subplot(2,1,2);qqplot(ptr_e22t);
title('QQ plot PTR AR-GARCH(2,2)-t Standardised Residuals');

% hist qq plot N-dist
figure;subplot(2,1,1);hist(ptr_e22tg,25);
title('Histogram of PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');
subplot(2,1,2);qqplot(ptr_e22tg);
title('QQ plot PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(ptr_e22tg, [10 15], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(ptr_e22tg.^2, [10 15], 0.05, [5 10])

%% 2. AR(p)-GARCH(2,2)-t
% Firstly consider AR(1)-GARCH(p,q) models with gaussian and student t errors
logL2=0;aict2=0;sict2=0;
for p=1:1
  for q=1:3
% NOTE that in MATLAB p is number of lagged conditional std dev's, q is number of lagged squared innovations
% this notation is different to lectures (p and q are switched) and most of the GARCH literature
%     Mdl = arima('ARLags',1,'Variance',garch(p,q));Mdl.Distribution='Gaussian';
%     [EstMdl,EstParamCov,logL,info] = estimate(Mdl,BHPr,'display','off'); 
%     aic(p,q)=-2*logL+2*(p+q+1);sic1(p,q)=-2*logL+log(length(BHPr))*(p+q+1);
    Mdl_2 = arima(p,0,q); Mdl_2.Variance=garch(2,2); Mdl_2.Distribution='T';
    [EstMdl,EstParamCov,logL2,info] = estimate(Mdl_2,ptr_rtn_is,'display','off');
    aict2(q)=-2*logL2+2*(6+q); sict2(q)=-2*logL2+log(length(ptr_rtn_is))*(6+q);
  end
end 

% AIC&SIC plot
figure;plot(aict2,'b+-');title('AIC & SIC for PTR ARMA-GARCH-t models');hold on;plot(sict2,'r+-');
legend('AIC','SIC');

[aict_min,aict_min_index]=min(aict);
[sict_min,sict_min_index]=min(sict);

% optimal AR-GARCH-t model fit
Mdl_2 = arima(2,0,2); Mdl_2.Variance=garch(2,2); Mdl_2.Distribution='T';
[EstMdl,EstParamCov,logL,info] = estimate(Mdl_2,ptr_rtn_is);
[ptr_a1122t,ptr_v1122t,logL] = infer(EstMdl,ptr_rtn_is); ptr_s1122t=sqrt(ptr_v1122t); %calculate innovations
ptr_e1122t=ptr_a1122t./ptr_s1122t; ptr_e1122t_sqr=ptr_e1122t.^2;
% G-error transformation
ptr_df=EstMdl.Distribution.DoF;
ptr_e1122tg=norminv(tcdf(sqrt(ptr_df)/sqrt(ptr_df-2)*ptr_e1122t,ptr_df));

% Plot the innovations, conditional standard deviations and log returns above one another
figure;subplot(3,1,1);plot(date_is(2:end,1),ptr_a1122t);
% xlim([date_is(2) date_is(end)]);
ylim([min(ptr_a1122t)-1 max(ptr_a1122t)+1]);                      % set range of y-axis
title('PTR AR-GARCH(2,2)-t Innovations');
subplot(3,1,2);plot(date_is(2:end,1),ptr_s1122t);
% xlim([date_is(2) date_is(end)]);
ylim([0 max(ptr_s1122t)+1]);                      % set range of y-axis
title('PTR AR-GARCH(2,2)-t Conditional Standard Deviations');
subplot(3,1,3);plot(date_is(2:end,1),ptr_rtn_is);
% xlim([date_is(2) date_is(end)]);
ylim([min(ptr_rtn_is)-1 max(ptr_rtn_is)+1]);                      % set range of y-axis
title('PTR Log returns');

% AR-GARCH(2,2)-T Residual Diagnostics
% Cstd & log rtn
figure;plot(date_is(2:end,1),ptr_s1122t,'b');hold on;plot(date_is(2:end,1),ptr_rtn_is,'r');legend('PTR AR-GARCH(2,2)-t Conditional Std','PTR Log Rtn');
title('Conditional Standard Deviations for PTR AR-GARCH(2,2)-t & Log Rtn ');

% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),ptr_e1122t,'b'); title('PTR Standardised Residuals');
figure; autocorr(ptr_e1122t,20); title('PTR Residual ACF Plot');
figure;autocorr(ptr_e1122t_sqr,20); title('PTR Squared Residual ACF Plot');

% hist qq plot t-dist
figure;subplot(2,1,1);hist(ptr_e1122t,25);
title('Histogram of PTR AR-GARCH(2,2)-t Standardised Residuals');
subplot(2,1,2);qqplot(ptr_e1122t);
title('QQ plot PTR AR-GARCH(2,2)-t Standardised Residuals');

% hist qq plot N-dist
figure;subplot(2,1,1);hist(ptr_e1122tg,25);
title('Histogram of PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');
subplot(2,1,2);qqplot(ptr_e1122tg);
title('QQ plot PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(ptr_e1122tg, [10 15], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(ptr_e1122tg.^2, [10 15], 0.05, [5 10])

% Final Model Selection
logL2=0;aict2=0;sict2=0;
for p=1:2
  for q=1:2
% NOTE that in MATLAB p is number of lagged conditional std dev's, q is number of lagged squared innovations
% this notation is different to lectures (p and q are switched) and most of the GARCH literature
%     Mdl = arima('ARLags',1,'Variance',garch(p,q));Mdl.Distribution='Gaussian';
%     [EstMdl,EstParamCov,logL,info] = estimate(Mdl,BHPr,'display','off'); 
%     aic(p,q)=-2*logL+2*(p+q+1);sic1(p,q)=-2*logL+log(length(BHPr))*(p+q+1);
    Mdl = arima(p,0,q); Mdl.Variance=garch(2,2); Mdl.Distribution='T';
    [EstMdl,EstParamCov,logL2,info] = estimate(Mdl,ptr_rtn_is,'display','off');
    aict2(p,q)=-2*logL2+2*(5+q+p); sict2(p,q)=-2*logL2+log(length(ptr_rtn_is))*(5+p+q);
  end
end

% AIC&SIC plot
figure;plot(aict2,'b+-');title('AIC & SIC for PTR ARMA-GARCH-t models');hold on;plot(sict2,'r+-');
legend('AIC','SIC');
% Trust AIC here ARMA(2,2)-GARCH(2,2)-T
%% PTR ARMA(2,2)-GARCH(2,2)-T
% optimal AR-GARCH-t model fit
Mdl = arima(2,0,2); Mdl.Variance=garch(2,2); Mdl.Distribution='T';
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,ptr_rtn_is);
[ptr_a2222t,ptr_v2222t,logL] = infer(EstMdl,ptr_rtn_is); ptr_s2222t=sqrt(ptr_v2222t); %calculate innovations
ptr_e2222t=ptr_a2222t./ptr_s2222t; ptr_e2222t_sqr=ptr_e2222t.^2;
% G-error transformation
ptr_df=EstMdl.Distribution.DoF;
ptr_e2222tg=norminv(tcdf(sqrt(ptr_df)/sqrt(ptr_df-2)*ptr_e2222t,ptr_df));




% Plot the innovations, conditional standard deviations and log returns above one another
figure;subplot(3,1,1);plot(date_is(2:end,1),ptr_a2222t);
% xlim([date_is(2) date_is(end)]);
ylim([min(ptr_a2222t)-1 max(ptr_a2222t)+1]);                      % set range of y-axis
title('PTR AR-GARCH(2,2)-t Innovations');
subplot(3,1,2);plot(date_is(2:end,1),ptr_s2222t);
% xlim([date_is(2) date_is(end)]);
ylim([0 max(ptr_s2222t)+1]);                      % set range of y-axis
title('PTR AR-GARCH(2,2)-t Conditional Standard Deviations');
subplot(3,1,3);plot(date_is(2:end,1),ptr_rtn_is);
% xlim([date_is(2) date_is(end)]);
ylim([min(ptr_rtn_is)-1 max(ptr_rtn_is)+1]);                      % set range of y-axis
title('PTR Log returns');

% ARMA(2,2)-GARCH(2,2)-T Residual Diagnostics
% Cstd & log rtn
figure;plot(date_is(2:end,1),ptr_s2222t,'b');hold on;plot(date_is(2:end,1),ptr_rtn_is,'r');legend('PTR AR-GARCH(2,2)-t Conditional Std','PTR Log Rtn');
title('Conditional Standard Deviations for PTR ARMA(2,2)-GARCH(2,2)-t & Log Rtn ');

% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),ptr_e2222t,'b'); title('PTR ARMA(2,2)-GARCH(2,2)-t Standardised Residuals');
figure; autocorr(ptr_e2222t,20); title('PTR ARMA(2,2)-GARCH(2,2)-t Residual ACF Plot');
figure;autocorr(ptr_e2222t_sqr,20); title('PTR ARMA(2,2)-GARCH(2,2)-t Squared Residual ACF Plot');

% hist qq plot t-dist
figure;subplot(2,1,1);hist(ptr_e1122t,25);
title('Histogram of PTR AR-GARCH(2,2)-t Standardised Residuals');
subplot(2,1,2);qqplot(ptr_e1122t);
title('QQ plot PTR AR-GARCH(2,2)-t Standardised Residuals');

% hist qq plot N-dist
figure;subplot(2,1,1);hist(ptr_e1122tg,25);
title('Histogram of PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');
subplot(2,1,2);qqplot(ptr_e1122tg);
title('QQ plot PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');

% JB test
[skewness(ptr_e1122tg) kurtosis(ptr_e1122tg)]
[h, p] = jbtest(ptr_e1122tg)

%% GJR-GARCH(2,2)-T
Mdl = gjr(2,2);Mdl.Distribution='t';Mdl.Offset=NaN; % specify model
[EstMdl,EstParamCov,LLF,info] = estimate(Mdl,ptr_rtn_is);
ptr_v22g = infer(EstMdl,ptr_rtn_is);ptr_s22g=sqrt(ptr_v22g); %infer the conditional variance and calculate standard deviations innovations
ptr_a22g = ptr_rtn_is-EstMdl.Offset; % innovations
ptr_e22g=ptr_a22g./ptr_s22g; % standardised residuals
% G-error transformation
ptr_df=EstMdl.Distribution.DoF;
ptr_e22gg=norminv(tcdf(sqrt(ptr_df)/sqrt(ptr_df-2)*ptr_e22g,ptr_df));

% Plot the innovations, conditional standard deviations and log returns above one another
figure;subplot(3,1,1);plot(date_is(2:end,1),ptr_e22g);
% xlim([date_is(2) date_is(end)]);
ylim([min(ptr_a22g)-1 max(ptr_a22g)+1]);                      % set range of y-axis
title('PTR GJR-GARCH(2,2)-t Innovations');
subplot(3,1,2);plot(date_is(2:end,1),ptr_s22g);
% xlim([date_is(2) date_is(end)]);
ylim([0 max(ptr_s22g)+1]);                      % set range of y-axis
title('PTR GJR-GARCH(2,2)-t Conditional Standard Deviations');
subplot(3,1,3);plot(date_is(2:end,1),ptr_rtn_is);
% xlim([date_is(2) date_is(end)]);
ylim([min(ptr_rtn_is)-1 max(ptr_rtn_is)+1]);                      % set range of y-axis
title('PTR Log returns');

% AR-GARCH(2,2)-T Residual Diagnostics
% Cstd & log rtn
figure;hold on;plot(date_is(2:end,1),ptr_rtn_is,'r');plot(date_is(2:end,1),ptr_s22g,'b');legend('PTR Log Rtn','PTR GJR-GARCH(2,2)-t Conditional Std');
title('Conditional Standard Deviations for PTR GJR-GARCH(2,2)-t & Log Rtn ');

% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),ptr_e22g,'b'); title('PTR GJR-GARCH(2,2)-t Standardised Residuals');
figure; autocorr(ptr_e22g,20); title('PTR GJR-GARCH(2,2)-t Residual ACF Plot');
figure;autocorr(ptr_e22g.^2,20); title('PTR GJR-GARCH(2,2)-t Squared Residual ACF Plot');

% hist qq plot t-dist
figure;subplot(2,1,1);hist(ptr_e22g,25);
title('Histogram of PTR AR-GARCH(2,2)-t Standardised Residuals');
subplot(2,1,2);qqplot(ptr_e22g);
title('QQ plot PTR AR-GARCH(2,2)-t Standardised Residuals');

% hist qq plot N-dist
figure;subplot(2,1,1);hist(ptr_e22gg,25);
title('Histogram of PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');
subplot(2,1,2);qqplot(ptr_e22gg);
title('QQ plot PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(ptr_e22gg, [10 15], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(ptr_e22gg.^2, [10 15], 0.05, [5 10])
%% ARMA(1,1)-GJR-GARCH(2,2)-T
Mdl = arima(1,0,1); Mdl.Variance=gjr(1,1); Mdl.Distribution='T';
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,ptr_rtn_is);

[ptr_a122g,ptr_v122g,logL] = infer(EstMdl,ptr_rtn_is); ptr_s122g=sqrt(ptr_v122g); %calculate innovations
ptr_e122g=ptr_a122g./ptr_s122g; ptr_e122g_sqr=ptr_e122g.^2;
% G-error transformation
ptr_df=EstMdl.Distribution.DoF;
ptr_e122gg=norminv(tcdf(sqrt(ptr_df)/sqrt(ptr_df-2)*ptr_e122g,ptr_df));



% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(ptr_e122gg, [10 15], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(ptr_e122gg.^2, [10 15], 0.05, [5 10])

logL2=0;aict2=0;sict2=0;
for p=1:2
  for q=1:2
% NOTE that in MATLAB p is number of lagged conditional std dev's, q is number of lagged squared innovations
% this notation is different to lectures (p and q are switched) and most of the GARCH literature
%     Mdl = arima('ARLags',1,'Variance',garch(p,q));Mdl.Distribution='Gaussian';
%     [EstMdl,EstParamCov,logL,info] = estimate(Mdl,BHPr,'display','off'); 
%     aic(p,q)=-2*logL+2*(p+q+1);sic1(p,q)=-2*logL+log(length(BHPr))*(p+q+1);
    Mdl = arima(1,0,1); Mdl.Variance=gjr(p,q); Mdl.Distribution='T';
    [EstMdl,EstParamCov,logL2,info] = estimate(Mdl,ptr_rtn_is,'display','off');
    aict2(p,q)=-2*logL2+2*(3+q+p); sict2(p,q)=-2*logL2+log(length(ptr_rtn_is))*(3+p+q);
  end
end

% AIC&SIC plot
figure;plot(aict2,'b+-');title('AIC & SIC for PTR ARMA-GARCH-t models');hold on;plot(sict2,'r+-');
legend('AIC','SIC');

%%  ARMA(2,2)-ARCH(10)-t
logL2=0;aict2=0;sict2=0;
for q=1:20
    Mdl = arima(2,0,2); Mdl.Variance=garch(0,q); Mdl.Distribution='T';
    [EstMdl,EstParamCov,logL2,info] = estimate(Mdl,ptr_rtn_is,'display','off');
    aict2(q)=-2*logL2+2*(5+q); sict2(q)=-2*logL2+log(length(ptr_rtn_is))*(5+q);
end

% AIC&SIC plot
figure;plot(aict2,'b+-');title('AIC & SIC for PTR ARMA-ARCH-t models');hold on;plot(sict2,'r+-');
legend('AIC','SIC');


Mdl = arima(2,0,2); Mdl.Variance=garch(0,10); Mdl.Distribution='T';
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,ptr_rtn_is);

[ptr_a122g,ptr_v122g,logL] = infer(EstMdl,ptr_rtn_is); ptr_s122g=sqrt(ptr_v122g); %calculate innovations
ptr_e122g=ptr_a122g./ptr_s122g; ptr_e122g_sqr=ptr_e122g.^2;
% G-error transformation
ptr_df=EstMdl.Distribution.DoF;
ptr_e122gg=norminv(tcdf(sqrt(ptr_df)/sqrt(ptr_df-2)*ptr_e122g,ptr_df));

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(ptr_e122gg, [10 15], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(ptr_e122gg.^2, [10 15], 0.05, [5 10])

% JB test
[skewness(ptr_e122gg) kurtosis(ptr_e122gg)]
[h, p] = jbtest(ptr_e122gg)

%% PTR GJR-GARCH-t 
for p=1:5
    for q=1:5
      Mdl = arima(1,0,0); Mdl.Variance=gjr(p,q); Mdl.Distribution='T';  
      [EstMdl,EstParamCov,logL,info] = estimate(Mdl,ptr_rtn_is);  
      aict(p,q)=-2*logL+2*(p+q+2); sict(p,q)=-2*logL+log(length(ptr_rtn_is))*(p+q+2);
    end
end

% AIC SIF Plot
figure;plot(1:5,aict(1,:),'b+-',1:5,aict(2,:),'k+-',1:5,aict(3,:),'r+-',1:5,aict(4,:),'g+-',1:5,aict(5,:),'y+-');
title('AIC & SIC for ptr AR(1)-GJR-GARCH-t Models');hold on;
plot(1:5,sict(1,:),'b+-',1:5,sict(2,:),'k+-',1:5,sict(3,:),'r+-',1:5,sict(4,:),'g+-',1:5,sict(5,:),'y+-');
legend('AIC,P=1','AIC,P=2','AIC,P=3','AIC,P=4','AIC,P=5','SIC,P=1','SIC,P=2','SIC,P=3','SIC,P=4','SIC,P=5');
xlabel('Number of Q');

% Optimal Model Training
Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T'; 
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,ptr_rtn_is);
[ptr_a111g,ptr_v111g,logL] = infer(EstMdl,ptr_rtn_is); ptr_s111g=sqrt(ptr_v111g); %calculate innovations
ptr_e111g=ptr_a111g./ptr_s111g; ptr_e111g_sqr=ptr_e111g.^2;
% G-error transformation
ptr_df=EstMdl.Distribution.DoF;
ptr_e111gg=norminv(tcdf(sqrt(ptr_df)/sqrt(ptr_df-2)*ptr_e111g,ptr_df));
% ptr_sigGJt=ptr_s11j; ptr_aGJt=ptr_a11j;

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(ptr_e111gg, [9 14], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(ptr_e111gg.^2, [9 14], 0.05, [5 10])

% JB Test
[skewness(ptr_e111gg) kurtosis(ptr_e111gg)]
[h, p] = jbtest(ptr_e111gg)

%% GJR-GARCH-t MMM
LLF=0;aict=0;sict=0;
for p=1:5
    for q = 1:5
        Mdl = gjr(p,q);Mdl.Distribution='t';Mdl.Offset=NaN; % specify model
        [EstMdl,EstParamCov,LLF,info] = estimate(Mdl,mmm_rtn_is);
        aict(p,q)=-2*LLF+2*(p+q+1); sict(p,q)=-2*LLF+log(length(mmm_rtn_is))*(p+q+1);
    end
end

figure;plot(aict,'b+-');title('AIC & SIC for PTR ARMA-ARCH-t models');hold on;plot(sict,'r+-');
legend('AIC','SIC');

Mdl = gjr(1,1);Mdl.Distribution='t';Mdl.Offset=NaN; % specify model
[EstMdl,EstParamCov,LLF,info] = estimate(Mdl,mmm_rtn_is);
mmm_v11j = infer(EstMdl,mmm_rtn_is);mmm_s11j=sqrt(mmm_v11j); %infer the conditional variance and calculate standard deviations innovations
mmm_a11j = mmm_rtn_is-EstMdl.Offset; % innovations
mmm_e11j=mmm_a11j./mmm_s11j; % standardised residuals
mmm_sigGJt=mmm_s11j; mmm_aGJt=mmm_a11j;

mmm_df=EstMdl.Distribution.DoF;
mmm_e11jg=norminv(tcdf(sqrt(mmm_df)/sqrt(mmm_df-2)*mmm_e11j,mmm_df));

% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),mmm_e11j,'b'); title('mmm ARMA(2,2)-GARCH(2,2)-t Standardised Residuals');
figure; autocorr(mmm_e11j,20); title('mmm ARMA(2,2)-GARCH(2,2)-t Residual ACF Plot');
figure;autocorr(mmm_e11j.^2,20); title('mmm ARMA(2,2)-GARCH(2,2)-t Squared Residual ACF Plot');


% NIC CURVE with a
% Get required estimated GJR model coefficients
mmm_a1GJt=cell2mat(EstMdl.ARCH);mmm_b1GJt=cell2mat(EstMdl.GARCH);
mmm_g1GJt=cell2mat(EstMdl.Leverage);
mmm_dfGJt=EstMdl.Distribution.DoF;  
mmm_a0GJt=EstMdl.Constant; 


%average volatility in each regime
[mmm_a0GJt/(1-mmm_a1GJt-mmm_g1GJt-mmm_b1GJt)  mmm_a0GJt/(1-mmm_a1GJt-mmm_b1GJt)]

% persistense
[mmm_a1GJt+mmm_g1GJt+mmm_b1GJt mmm_a1GJt+mmm_b1GJt]

% Caculate values for NIC plot
mmm_a=min(mmm_aGJt):0.01:max(mmm_aGJt);   % range of plot
mmm_sg=var(mmm_aGJt);	% sample variance of innovations
mmm_sigt=mmm_a0GJt+mmm_a1GJt*mmm_a.^2+mmm_b1GJt*mmm_sg+(mmm_g1GJt.*(mmm_a<0)).*(mmm_a.^2); % asymmetric curve values
mmm_sigt2=mmm_a0GJt+(mmm_a1GJt+mmm_g1GJt/2)*mmm_a.^2+mmm_b1GJt*mmm_sg;	 % symmetric curve values

% NIC plot for a(t-1)
figure;plot(mmm_a,mmm_sigt);axis([min(mmm_aGJt)-0.5 max(mmm_aGJt)+0.5 0 max(mmm_sigt)+1]);
hold on;plot(mmm_a,mmm_sigt2,'r--');
title('ARMA-GJR-t NIC curve for a(t-1)');


% compare precited conditional standard deviation values following positive
% and negative shocks of size 2
mmm_a=-1.746633160921687;mmm_sigtm2=mmm_a0GJt+mmm_a1GJt*mmm_a^2+mmm_b1GJt*mmm_sg+(mmm_g1GJt*(mmm_a<0))*(mmm_a^2); 
mmm_a=1.746633160921687;mmm_sigta2=mmm_a1GJt+mmm_a1GJt*mmm_a^2+mmm_b1GJt*mmm_sg+(mmm_g1GJt*(mmm_a<0))*(mmm_a^2);
[mmm_sigtm2 mmm_sigta2 mmm_sigtm2/mmm_sigta2] 


% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(mmm_e11jg, [10 15], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(mmm_e11jg.^2, [10 15], 0.05, [5 10])


% hist qq plot t-dist
figure;subplot(2,1,1);hist(mmm_e11j,25);
title('Histogram of PTR AR-GARCH(2,2)-t Standardised Residuals');
subplot(2,1,2);qqplot(mmm_e11j);
title('QQ plot PTR AR-GARCH(2,2)-t Standardised Residuals');

% hist qq plot N-dist
figure;subplot(2,1,1);hist(mmm_e11jg,25);
title('Histogram of PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');
subplot(2,1,2);qqplot(mmm_e11jg);
title('QQ plot PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');

% JB test
[skewness(mmm_e11jg) kurtosis(mmm_e11jg)]
[h, p] = jbtest(mmm_e11jg)

var(mmm_rtn_is)

%% ARMA-GJR-GARCH-t MMM
Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T';
[EstMdl,EstParamCov,logL,info] = estimate(Mdl,mmm_rtn_is);

[mmm_a31g,mmm_v31g,logL] = infer(EstMdl,mmm_rtn_is); mmm_s31g=sqrt(mmm_v31g); %calculate innovations
mmm_e31g=mmm_a31g./mmm_s31g; mmm_e31g_sqr=mmm_e31g.^2;
% G-error transformation
mmm_df=EstMdl.Distribution.DoF;
mmm_e31gg=norminv(tcdf(sqrt(mmm_df)/sqrt(mmm_df-2)*mmm_e31g,mmm_df));


% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),mmm_e31g,'b'); title('mmm ARMA(2,2)-GARCH(2,2)-t Standardised Residuals');
figure; autocorr(mmm_e31g,20); title('mmm ARMA(2,2)-GARCH(2,2)-t Residual ACF Plot');
figure;autocorr(mmm_e31g.^2,20); title('mmm ARMA(2,2)-GARCH(2,2)-t Squared Residual ACF Plot');

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(mmm_e31gg, [10 15], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(mmm_e31gg.^2, [10 15], 0.05, [5 10])


%% EGARCH(3,1)
LLF=0;aict=0;sict=0;
for p=1:4
    for q = 1:4
        Mdl = egarch(p-1,q);Mdl.Distribution='t';Mdl.Offset=NaN; % specify model
        [EstMdl,EstParamCov,LLF,info] = estimate(Mdl,mmm_rtn_is);
        aict(p,q)=-2*LLF+2*(p+q); sict(p,q)=-2*LLF+log(length(mmm_rtn_is))*(p+q);
    end
end

figure;plot(aict,'b+-');title('AIC & SIC for PTR ARMA-ARCH-t models');hold on;plot(sict,'r+-');
legend('AIC','SIC');


Mdl = egarch(3,1);Mdl.Distribution='t';Mdl.Offset=NaN; % specify model
[EstMdl,EstParamCov,LLF,info] = estimate(Mdl,mmm_rtn_is);
mmm_v31e = infer(EstMdl,mmm_rtn_is);mmm_s31e=sqrt(mmm_v31e); %infer the conditional variance and calculate standard deviations innovations
mmm_a31e = mmm_rtn_is-EstMdl.Offset; % innovations
mmm_e31e=mmm_a31e./mmm_s31e; % standardised residuals

% G-error transformation
mmm_df=EstMdl.Distribution.DoF;
mmm_e31eg=norminv(tcdf(sqrt(mmm_df)/sqrt(mmm_df-2)*mmm_e31e,mmm_df));


% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),mmm_e31e,'b'); title('mmm ARMA(2,2)-GARCH(2,2)-t Standardised Residuals');
figure; autocorr(mmm_e31e,20); title('mmm ARMA(2,2)-GARCH(2,2)-t Residual ACF Plot');
figure;autocorr(mmm_e31e.^2,20); title('mmm ARMA(2,2)-GARCH(2,2)-t Squared Residual ACF Plot');


% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(mmm_e31eg, [10 15], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(mmm_e31eg.^2, [10 15], 0.05, [5 10])

% hist qq plot t-dist
figure;subplot(2,1,1);hist(mmm_e31e,25);
title('Histogram of PTR AR-GARCH(2,2)-t Standardised Residuals');
subplot(2,1,2);qqplot(mmm_e31e);
title('QQ plot PTR AR-GARCH(2,2)-t Standardised Residuals');

% hist qq plot N-dist
figure;subplot(2,1,1);hist(mmm_e31eg,25);
title('Histogram of PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');
subplot(2,1,2);qqplot(mmm_e31eg);
title('QQ plot PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');

% JB test
[skewness(mmm_e31eg) kurtosis(mmm_e31eg)]
[h, p] = jbtest(mmm_e31eg)

%% ARMA-EGARCH(3,1)-t

for p=1:4
    for q = 1:4
        Mdl = egarch(p-1,q);Mdl.Distribution='t';Mdl.Offset=NaN; % specify model
        [EstMdl,EstParamCov,LLF,info] = estimate(Mdl,mmm_rtn_is);
        aict(p,q)=-2*LLF+2*(p+q); sict(p,q)=-2*LLF+log(length(mmm_rtn_is))*(p+q);
    end
end

Mdleg = arima('ARLags',1,'Variance',egarch(3,1));
Mdleg.Distribution='t'; % specify model
[EstMdl,EstParamCov,LLF,info] = estimate(Mdleg,mmm_rtn_is);  % Estimate model
[mmm_aeg31t,mmm_veg31t,LLF] = infer(EstMdl,mmm_rtn_is);mmm_seg31t=sqrt(mmm_veg31t); %infer the conditional variance and calculate standard deviations innovations
mmm_eeg31t=mmm_aeg31t./mmm_seg31t;
% G-error transformation
mmm_df=EstMdl.Distribution.DoF;
mmm_eeg31tg=norminv(tcdf(sqrt(mmm_df)/sqrt(mmm_df-2)*mmm_eeg31t,mmm_df));

% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),mmm_eeg31t,'b'); title('mmm ARMA(2,2)-GARCH(2,2)-t Standardised Residuals');
figure; autocorr(mmm_eeg31t,20); title('mmm ARMA(2,2)-GARCH(2,2)-t Residual ACF Plot');
figure;autocorr(mmm_eeg31t.^2,20); title('mmm ARMA(2,2)-GARCH(2,2)-t Squared Residual ACF Plot');

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(mmm_eeg31tg, [10 15], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(mmm_eeg31tg.^2, [10 15], 0.05, [5 10])

%% MMM - IGARCH
[PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS]=igarch(mmm_rtn_is,1,1,[],[],0);
SHT=sqrt(HT);
mmm_eHT=mmm_rtn_is./SHT;

% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),mmm_eHT,'b'); title('mmm ARMA(2,2)-GARCH(2,2)-t Standardised Residuals');
figure; autocorr(mmm_eHT,20); title('mmm ARMA(2,2)-GARCH(2,2)-t Residual ACF Plot');
figure;autocorr(mmm_eHT.^2,20); title('mmm ARMA(2,2)-GARCH(2,2)-t Squared Residual ACF Plot');


% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(mmm_eHT, [7 12], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(mmm_eHT.^2, [7 12], 0.05, [5 10])

% hist qq plot t-dist
figure;subplot(2,1,1);hist(mmm_eHT,25);
title('Histogram of PTR AR-GARCH(2,2)-t Standardised Residuals');
subplot(2,1,2);qqplot(mmm_eHT);
title('QQ plot PTR AR-GARCH(2,2)-t Standardised Residuals');

% hist qq plot N-dist
figure;subplot(2,1,1);hist(mmm_eHT,25);
title('Histogram of PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');
subplot(2,1,2);qqplot(mmm_eHT);
title('QQ plot PTR AR-GARCH(2,2)-t Normal Inversed Standardised Residuals');

% JB test
[skewness(mmm_eHT) kurtosis(mmm_eHT)]
[h, p] = jbtest(mmm_eHT)
%% Coke GJR-GARCH(1,1)-t
Mdl = gjr(1,1);Mdl.Distribution='t';Mdl.Offset=NaN; % specify model
[EstMdl,EstParamCov,LLF,info] = estimate(Mdl,coke_rtn_is);
coke_v11j = infer(EstMdl,coke_rtn_is);coke_s11j=sqrt(coke_v11j); %infer the conditional variance and calculate standard deviations innovations
coke_a11j = coke_rtn_is-EstMdl.Offset; % innovations
coke_e11j=coke_a11j./coke_s11j; % standardised residuals
coke_sigGJt=coke_s11j; coke_aGJt=coke_a11j;

coke_df=EstMdl.Distribution.DoF;
coke_e11jg=norminv(tcdf(sqrt(coke_df)/sqrt(coke_df-2)*coke_e11j,coke_df));

% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),coke_e11j,'b'); title('coke ARMA(2,2)-GARCH(2,2)-t Standardised Residuals');
figure; autocorr(coke_e11j,20); title('coke ARMA(2,2)-GARCH(2,2)-t Residual ACF Plot');
figure;autocorr(coke_e11j.^2,20); title('coke ARMA(2,2)-GARCH(2,2)-t Squared Residual ACF Plot');

% NIC CURVE with a
% Get required estimated GJR model coefficients
coke_a1GJt=cell2mat(EstMdl.ARCH);coke_b1GJt=cell2mat(EstMdl.GARCH);
coke_g1GJt=cell2mat(EstMdl.Leverage);
coke_dfGJt=EstMdl.Distribution.DoF;  
coke_a0GJt=EstMdl.Constant; 


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
title('ARMA-GJR-t NIC curve for a(t-1)');


% compare precited conditional standard deviation values following positive
% and negative shocks of size 2
coke_a=-2;coke_sigtm2=coke_a0GJt+coke_a1GJt*coke_a^2+coke_b1GJt*coke_sg+(coke_g1GJt*(coke_a<0))*(coke_a^2); 
coke_a=2;coke_sigta2=coke_a1GJt+coke_a1GJt*coke_a^2+coke_b1GJt*coke_sg+(coke_g1GJt*(coke_a<0))*(coke_a^2);
[coke_sigtm2 coke_sigta2 coke_sigtm2/coke_sigta2] 
%% Coke AR(1)-GJR-GARCH(1,1)-t
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
coke_sigGJt=coke_s11j; coke_aGJt=coke_a11j;

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
title('Conditional Std for Coke AR(1)-GJR-GARCH(1,1)-t & Log Rtn ');

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

%% Coke-GARCH-t
logL=0;aict=0;sict=0;
for p=1:5
    for q=1:5
      Mdl = arima(1,0,0); Mdl.Variance=garch(p,q); Mdl.Distribution='T';  
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

% AIC SIF Plot
figure;plot(1:5,aict(1,:),'b+-',1:5,aict(2,:),'k+-',1:5,aict(3,:),'r+-',1:5,aict(4,:),'g+-',1:5,aict(5,:),'y+-');
title('AIC & SIC for Coke AR(1)-GJR-GARCH-t Models');hold on;
plot(1:5,sict(1,:),'b+-',1:5,sict(2,:),'k+-',1:5,sict(3,:),'r+-',1:5,sict(4,:),'g+-',1:5,sict(5,:),'y+-');
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
coke_sigGJt=coke_s11j; coke_aGJt=coke_a11j;

% Plot the innovations, conditional standard deviations and log returns above one another
figure;subplot(3,1,1);plot(date_is(2:end,1),coke_a111t);
% xlim([date_is(2) date_is(end)]);
ylim([min(coke_a111t)-1 max(coke_a111t)+1]);                      % set range of y-axis
title('Coke AR(1)-GARCH(1,1)-t Innovations');
subplot(3,1,2);plot(date_is(2:end,1),coke_s111t);
% xlim([date_is(2) date_is(end)]);
ylim([0 max(coke_s111t)+1]);                      % set range of y-axis
title('Coke AR(1)-GARCH(1,1)-t Conditional Standard Deviations');
subplot(3,1,3);plot(date_is(2:end,1),coke_rtn_is);
% xlim([date_is(2) date_is(end)]);
ylim([min(coke_rtn_is)-1 max(coke_rtn_is)+1]);                      % set range of y-axis
title('Coke Log returns');

% Shoulder Plot
figure;plot(date_is(2:end,1),coke_s111t,'b');hold on;plot(date_is(2:end,1),coke_rtn_is,'r');legend('Coke AR(1)-GARCH(1,1)-t Conditional Std','Coke Log Rtn');
title('Conditional Std for Coke AR(1)-GARCH(1,1)-t & Log Rtn ');


% ACF ACF Squared plot
% Standardised Resid Plot
figure;subplot(3,1,1); plot(date_is(2:end,1),coke_e111t,'b'); title('Coke AR(1)-GARCH(1,1)-t Standardised Residuals');
figure; autocorr(coke_e111t,20); title('Coke AR(1)-GARCH(1,1)-t Residual ACF Plot');
figure;autocorr(coke_e111t_sqr,20); title('Coke AR(1)-GARCH(1,1)-t Squared Residual ACF Plot');


% hist qq plot t-dist
figure;subplot(2,1,1);hist(coke_e111t,25);
title('AR(1)-GARCH(1,1)-t Standardised Residuals');
subplot(2,1,2);qqplot(coke_e111t);
title('QQ plot Coke AR(1)-GARCH(1,1)-t Standardised Residuals');

% hist qq plot N-dist
figure;subplot(2,1,1);hist(coke_e111tg,25);
title('Coke AR(1)-GARCH(1,1)-t Normal Inversed Standardised Residuals');
subplot(2,1,2);qqplot(coke_e111tg);
title('Coke AR(1)-GARCH(1,1)-t Normal Inversed Standardised Residuals');

% LB test on transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(coke_e111tg, [9 14], 0.05, [5 10])
%LB test on squared transformed standardised residuals
[H, pValue, Qstat, CriticalValue] = lbqtest(coke_e111tg.^2, [9 14], 0.05, [5 10])

% JB Test
[skewness(coke_e111tg) kurtosis(coke_e111tg)]
[h, p] = jbtest(coke_e111tg)

%% I-GARCH(1,1)-N
[PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS]=igarch(coke_rtn_is,1,1,[],[],0);
SHT=sqrt(HT);
coke_eHT=coke_rtn_is./SHT;

% Shoulder Plot
figure;plot(date_is(2:end,1),SHT,'b');hold on;plot(date_is(2:end,1),coke_rtn_is,'r');legend('Coke AR(1)-GJR-GARCH(1,1)-t Conditional Std','Coke Log Rtn');
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
