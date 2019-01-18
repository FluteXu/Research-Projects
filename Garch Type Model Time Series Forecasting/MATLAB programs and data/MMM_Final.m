%% 1. GJR-GARCH-t MMM
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
%% 2. MMM - IGARCH
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
%% MMM - Risk Metrics
 ht=igarch_core(mmm_rtn_is,[0 0.06 0.94],[],1,1,[],[],[],0);
