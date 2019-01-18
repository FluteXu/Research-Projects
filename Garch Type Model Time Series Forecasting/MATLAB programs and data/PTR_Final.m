%% 1. PTR ARMA(2,2)-GARCH(2,2)-T
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

%% 2.  ARMA(2,2)-ARCH(10)-t
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


