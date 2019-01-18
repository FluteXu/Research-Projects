
%% 1. Moving Horizon 1 forecast - Out of Sample
% coke
coke_return_forecast_1=zeros(length(coke_rtn_os),1); 
coke_return_forecast_2=zeros(length(coke_rtn_os),1); 
% HS25 - Return
window = 25;
weights = repelem(1/window, window);
coke_return_forecast_4 = tsmovavg(coke_rtn,'w',weights,1);
coke_return_forecast_4 =coke_return_forecast_4 (end-839:end);
% HS5 - Return
window = 5;
weights = repelem(1/window, window);
coke_return_forecast_5 = tsmovavg(coke_rtn,'w',weights,1);
coke_return_forecast_5=coke_return_forecast_5(end-839:end);

coke_vol_forecast_1=zeros(length(coke_rtn_os),1);
coke_vol_forecast_2=zeros(length(coke_rtn_os),1);
coke_vol_forecast_3= zeros(length(coke_rtn_os),1);
% HS25 Variance
coke_vol_forecast_4=coke_hs25(end-839:end);
% HS5 Variance
coke_vol_forecast_5=coke_hs5(end-839:end);                                                     


% msft
msft_return_forecast_1=zeros(length(msft_rtn_os),1); 
msft_return_forecast_2=zeros(length(msft_rtn_os),1); 
% HS25 - Return
window = 25;
weights = repelem(1/window, window);
msft_return_forecast_4 = tsmovavg(msft_rtn,'w',weights,1);
msft_return_forecast_4 =msft_return_forecast_4 (end-839:end);
% HS5 - Return
window = 5;
weights = repelem(1/window, window);
msft_return_forecast_5 = tsmovavg(msft_rtn,'w',weights,1);
msft_return_forecast_5=msft_return_forecast_5(end-839:end);

msft_vol_forecast_1=zeros(length(msft_rtn_os),1);
msft_vol_forecast_2=zeros(length(msft_rtn_os),1);
msft_vol_forecast_3= zeros(length(msft_rtn_os),1);
% HS25 Variance
msft_vol_forecast_4=msft_hs25(end-839:end);
% HS5 Variance
msft_vol_forecast_5=msft_hs5(end-839:end); 

% jnj
jnj_return_forecast_1=zeros(length(jnj_rtn_os),1); 
jnj_return_forecast_2=zeros(length(jnj_rtn_os),1); 
% HS25 - Return
window = 25;
weights = repelem(1/window, window);
jnj_return_forecast_4 = tsmovavg(jnj_rtn,'w',weights,1);
jnj_return_forecast_4 =jnj_return_forecast_4 (end-839:end);
% HS5 - Return
window = 5;
weights = repelem(1/window, window);
jnj_return_forecast_5 = tsmovavg(jnj_rtn,'w',weights,1);
jnj_return_forecast_5=jnj_return_forecast_5(end-839:end);

jnj_vol_forecast_1=zeros(length(jnj_rtn_os),1);
jnj_vol_forecast_2=zeros(length(jnj_rtn_os),1);
jnj_vol_forecast_3= zeros(length(jnj_rtn_os),1);
% HS25 Variance
jnj_vol_forecast_4=jnj_hs25(end-839:end);
% HS5 Variance
jnj_vol_forecast_5=jnj_hs5(end-839:end); 

% mmm
mmm_return_forecast_1=zeros(length(mmm_rtn_os),1); 
% HS25 - Return
window = 100;
weights = repelem(1/window, window);
mmm_return_forecast_3 = tsmovavg(mmm_rtn,'w',weights,1);
mmm_return_forecast_3 =mmm_return_forecast_3 (end-839:end);
% HS25 - Return
window = 25;
weights = repelem(1/window, window);
mmm_return_forecast_4 = tsmovavg(mmm_rtn,'w',weights,1);
mmm_return_forecast_4 =mmm_return_forecast_4 (end-839:end);
% HS5 - Return
window = 5;
weights = repelem(1/window, window);
mmm_return_forecast_5 = tsmovavg(mmm_rtn,'w',weights,1);
mmm_return_forecast_5=mmm_return_forecast_5(end-839:end);

mmm_vol_forecast_1=zeros(length(mmm_rtn_os),1);
mmm_vol_forecast_2=zeros(length(mmm_rtn_os),1);
% HS100 Variance
mmm_vol_forecast_3=mmm_hs100(end-839:end);
% HS25 Variance
mmm_vol_forecast_4=mmm_hs25(end-839:end);
% HS5 Variance
mmm_vol_forecast_5=mmm_hs5(end-839:end); 

% ptr
ptr_return_forecast_1=zeros(length(ptr_rtn_os),1); 
ptr_return_forecast_2=zeros(length(ptr_rtn_os),1); 
% HS25 - Return
window = 25;
weights = repelem(1/window, window);
ptr_return_forecast_4 = tsmovavg(ptr_rtn,'w',weights,1);
ptr_return_forecast_4 =ptr_return_forecast_4 (end-839:end);
% HS5 - Return
window = 5;
weights = repelem(1/window, window);
ptr_return_forecast_5 = tsmovavg(ptr_rtn,'w',weights,1);
ptr_return_forecast_5=ptr_return_forecast_5(end-839:end);

ptr_vol_forecast_1=zeros(length(ptr_rtn_os),1);
ptr_vol_forecast_2=zeros(length(ptr_rtn_os),1);
ptr_vol_forecast_3= zeros(length(ptr_rtn_os),1);
% HS25 Variance
ptr_vol_forecast_4=ptr_hs25(end-839:end);
% HS5 Variance
ptr_vol_forecast_5=ptr_hs5(end-839:end); 

for i = 1:length(coke_rtn_os)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Coke %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % return sample window
    coke_return_sample = coke_rtn(1:end-841+i);
    
    %%%%%%%%%%%%%%%%%%%%%%% 1. Coke AR(1)-GARCH(1,1)-t
    Mdl = arima(1,0,0); Mdl.Variance=garch(1,1); Mdl.Distribution='T';  % refit the oject everytime
    [EstMdl,EstParamCov,logL,info] = estimate(Mdl,coke_return_sample,'display','off'); 
    %%% forecast return & volatility
    [coke_return_forecast_1(i),~,coke_vol_forecast_1(i)]=forecast(EstMdl,1,'Y0',coke_return_sample);
    
    %%%%%%%%%%%%%%%%%%%% 2. Coke AR(1)-GJR-GARCH(1,1)-t
    Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T'; 
    [EstMdl,EstParamCov,logL,info] = estimate(Mdl,coke_return_sample,'display','off');
    %%% forecast return & volatility
    [coke_return_forecast_2(i),~,coke_vol_forecast_2(i)]=forecast(EstMdl,1,'Y0',coke_return_sample);

    %%%%%%%%%%%%%%%%%%%%%%%%%%% 3. I-GARCH(1,1)-N (For volatility only)
    [PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS]=igarch(coke_return_sample,1,1,[],[],0);
    %%% forecast volatility
    coke_vol_forecast_3(i)=PARAMETERS*coke_return_sample(end)^2+(1-PARAMETERS)*HT(end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% MSFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    msft_return_sample = msft_rtn(1:end-841+i);
    %%%%%%%%%%%%%%%%%%%%%%%%% 1. MSFT AR(1)-GARCH(1,1)-t
    Mdl = arima(1,0,0); Mdl.Variance=garch(1,1); Mdl.Distribution='T'; 
    [EstMdl,EstParamCov,logL,info] = estimate(Mdl,msft_return_sample,'display','off');

    %%% forecast return volatility
    [msft_return_forecast_1(i),~,msft_vol_forecast_1(i)]=forecast(EstMdl,1,'Y0',msft_return_sample);
    
    %%%%%%%%%%%%% 2. MSFT AR(1)-GJR-GARCH(1,1)-t
    Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T'; 
    [EstMdl,EstParamCov,logL,info] = estimate(Mdl,msft_return_sample,'display','off');
    %%% forecast return volatility
    [msft_return_forecast_2(i),~,msft_vol_forecast_2(i)]=forecast(EstMdl,1,'Y0',msft_return_sample);

    %%%%%%%%%%%%%% 3. MSFT-IGARCH-N (for risk only)
    [PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS]=igarch(msft_return_sample,1,1,[],[],0);
    
    %%% forecast volatility
    msft_vol_forecast_3(i)=PARAMETERS*msft_return_sample(end)^2+(1-PARAMETERS)*HT(end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% JNJ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    jnj_return_sample = jnj_rtn(1:end-841+i);
    %%%%%%%%%%%%%%%%%1.  JNP AR-GARCH-t
    Mdl = arima(1,0,0); Mdl.Variance=garch(1,1); Mdl.Distribution='T'; 
    [EstMdl,EstParamCov,logL,info] = estimate(Mdl,jnj_return_sample,'display','off');
    %%% forecast return volatility
    [jnj_return_forecast_1(i),~,jnj_vol_forecast_1(i)]=forecast(EstMdl,1,'Y0',jnj_return_sample);

    %%%%%%%%%%%%%%%%%2. JNP AR-GJR-GARCH-t
    Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T'; 
    [EstMdl,EstParamCov,logL,info] = estimate(Mdl,jnj_return_sample,'display','off');
    %%% forecast return volatility
    [jnj_return_forecast_2(i),~,jnj_vol_forecast_2(i)]=forecast(EstMdl,1,'Y0',jnj_return_sample);
    
    %%%%%%%%%%%%%%%%%3. JNJ-IGARCH-N (for risk only)
    [PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS]=igarch(jnj_return_sample,1,1,[],[],0);
    %%% forecast volatility
    jnj_vol_forecast_3(i)=PARAMETERS*jnj_return_sample(end)^2+(1-PARAMETERS)*HT(end);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% MMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mmm_return_sample = mmm_rtn(1:end-841+i);
    %%%%%%%%%%%%%%%%%%%% 1. MMM-GJR-GARCH-t MMM
     Mdl = arima(1,0,0); Mdl.Variance=garch(1,1); Mdl.Distribution='T'; 
    [EstMdl,EstParamCov,logL,info] = estimate(Mdl,mmm_return_sample,'display','off');
    % forecast return & volatility
    [mmm_return_forecast_1(i),~,mmm_vol_forecast_1(i)]=forecast(EstMdl,1,'Y0',mmm_return_sample);
    %%%%%%%%%%%%%%%%%%%% 2. MMM - IGARCH
    [PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS]=igarch(mmm_return_sample,1,1,[],[],0);
    %%% forecast volatility
    mmm_vol_forecast_2(i)=PARAMETERS*mmm_return_sample(end)^2+(1-PARAMETERS)*HT(end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% PTR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ptr_return_sample = ptr_rtn(1:end-841+i);
    %%%%%%%%%%%%%%%%%%%%%% 1. PTR ARMA(2,2)-GARCH(2,2)-T
    Mdl = arima(2,0,2); Mdl.Variance=garch(2,2); Mdl.Distribution='T';
    [EstMdl,EstParamCov,logL,info] = estimate(Mdl,ptr_return_sample,'display','off');
    %%% forecast return volatility
    [ptr_return_forecast_1(i),~,ptr_vol_forecast_1(i)]=forecast(EstMdl,1,'Y0',ptr_return_sample);
    
    %%%%%%%%%%%%%%%%%%%%%% 2. PTR ARMA(2,2)-ARCH(10)-t
    Mdl = arima(2,0,2); Mdl.Variance=garch(0,10); Mdl.Distribution='T';
    [EstMdl,EstParamCov,logL,info] = estimate(Mdl,ptr_return_sample,'display','off');
    %%% forecast return volatility
    [ptr_return_forecast_2(i),~,ptr_vol_forecast_2(i)]=forecast(EstMdl,1,'Y0',ptr_return_sample);
    
    %%%%%%%%%%%%%%%%%%%%%%%3. PTR IGARCH(1,1)-N (for risk only)
    [PARAMETERS,LL,HT,VCVROBUST,VCV,SCORES,DIAGNOSTICS]=igarch(ptr_return_sample,1,1,[],[],0);
    %%% forecast volatility
    ptr_vol_forecast_3(i)=PARAMETERS*ptr_return_sample(end)^2+(1-PARAMETERS)*HT(end);
    
end

%% 2. Accuracy Assessment
%%%%%%%%%%%%%%%%%%%% Proxy Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=840
nis=length(coke_rtn)-fs
%% coke proxy
% proxy 1
coke_prox1C=abs(coke_rtn-mean(coke_rtn));
%proxy 2
coke_highC=COKE(:,3);coke_lowC=COKE(:,4);coke_openC=COKE(:,2);coke_closeC=COKE(:,5);
coke_rangC=100*log(coke_highC./coke_lowC);coke_rangC(coke_rangC <= 0)=mean([0 min(coke_rangC(coke_rangC>0))]);
coke_prox2C=sqrt(0.3607*(coke_rangC.^2));coke_prox2C=coke_prox2C(2:end);
%proxy 3
coke_prox3C=1.107*coke_prox2C.^2 + 0.68*((100*log(coke_openC(2:end)./coke_closeC(1:end-1))).^2);
coke_prox3C=sqrt(coke_prox3C);
%proxy 4
coke_prox4C=exp(2*log(coke_rangC)-0.86+2*0.29^2);
coke_prox4C=sqrt(coke_prox4C(2:end));

%in sample and forecast proxy
coke_prox1Ci=coke_prox1C(1:nis);coke_prox2Ci=coke_prox2C(1:nis);
coke_prox3Ci=coke_prox3C(1:nis);coke_prox4Ci=coke_prox4C(1:nis);
coke_prox1Cf=coke_prox1C(nis+1:end); coke_prox2Cf=coke_prox2C(nis+1:end); 
coke_prox3Cf=coke_prox3C(nis+1:end); coke_prox4Cf=coke_prox4C(nis+1:end); 

% figure;
% plot(coke_prox1C);hold on;
% plot(coke_prox2C);plot(coke_prox3C);plot(coke_prox3C);

%% mmm proxy
% proxy 1
mmm_prox1C=abs(mmm_rtn-mean(mmm_rtn));
%proxy 2
mmm_highC=MMM(:,3);mmm_lowC=MMM(:,4);mmm_openC=MMM(:,2);mmm_closeC=MMM(:,5);
mmm_rangC=100*log(mmm_highC./mmm_lowC);mmm_rangC(mmm_rangC <= 0)=mean([0 min(mmm_rangC(mmm_rangC>0))]);
mmm_prox2C=sqrt(0.3607*(mmm_rangC.^2));mmm_prox2C=mmm_prox2C(2:end);
%proxy 3
mmm_prox3C=1.107*mmm_prox2C.^2 + 0.68*((100*log(mmm_openC(2:end)./mmm_closeC(1:end-1))).^2);
mmm_prox3C=sqrt(mmm_prox3C);
%proxy 4
mmm_prox4C=exp(2*log(mmm_rangC)-0.86+2*0.29^2);
mmm_prox4C=sqrt(mmm_prox4C(2:end));


%in sample and forecast proxy
mmm_prox1Ci=mmm_prox1C(1:nis);mmm_prox2Ci=mmm_prox2C(1:nis);
mmm_prox3Ci=mmm_prox3C(1:nis);mmm_prox4Ci=mmm_prox4C(1:nis);
mmm_prox1Cf=mmm_prox1C(nis+1:end); mmm_prox2Cf=mmm_prox2C(nis+1:end); 
mmm_prox3Cf=mmm_prox3C(nis+1:end); mmm_prox4Cf=mmm_prox4C(nis+1:end); 

% figure;
% plot(mmm_prox1C);hold on;
% plot(mmm_prox2C);plot(mmm_prox3C);plot(mmm_prox3C);

%% jnj proxy
% proxy 1
jnj_prox1C=abs(jnj_rtn-mean(jnj_rtn));
%proxy 2
jnj_highC=JNJ(:,3);jnj_lowC=JNJ(:,4);jnj_openC=JNJ(:,2);jnj_closeC=JNJ(:,5);
jnj_rangC=100*log(jnj_highC./jnj_lowC);jnj_rangC(jnj_rangC <= 0)=mean([0 min(jnj_rangC(jnj_rangC>0))]);
jnj_prox2C=sqrt(0.3607*(jnj_rangC.^2));jnj_prox2C=jnj_prox2C(2:end);
%proxy 3
jnj_prox3C=1.107*jnj_prox2C.^2 + 0.68*((100*log(jnj_openC(2:end)./jnj_closeC(1:end-1))).^2);
jnj_prox3C=sqrt(jnj_prox3C);
%proxy 4
jnj_prox4C=exp(2*log(jnj_rangC)-0.86+2*0.29^2);
jnj_prox4C=sqrt(jnj_prox4C(2:end));


%in sample and forecast proxy
jnj_prox1Ci=jnj_prox1C(1:nis);jnj_prox2Ci=jnj_prox2C(1:nis);
jnj_prox3Ci=jnj_prox3C(1:nis);jnj_prox4Ci=jnj_prox4C(1:nis);
jnj_prox1Cf=jnj_prox1C(nis+1:end); jnj_prox2Cf=jnj_prox2C(nis+1:end); 
jnj_prox3Cf=jnj_prox3C(nis+1:end); jnj_prox4Cf=jnj_prox4C(nis+1:end); 

% figure;
% plot(jnj_prox1C);hold on;
% plot(jnj_prox2C);plot(jnj_prox3C);plot(jnj_prox3C);

%% msft proxy
% proxy 1
msft_prox1C=abs(msft_rtn-mean(msft_rtn));
%proxy 2
msft_highC=MSFT(:,3);msft_lowC=MSFT(:,4);msft_openC=MSFT(:,2);msft_closeC=MSFT(:,5);
msft_rangC=100*log(msft_highC./msft_lowC);msft_rangC(msft_rangC <= 0)=mean([0 min(msft_rangC(msft_rangC>0))]);
msft_prox2C=sqrt(0.3607*(msft_rangC.^2));msft_prox2C=msft_prox2C(2:end);
%proxy 3
msft_prox3C=1.107*msft_prox2C.^2 + 0.68*((100*log(msft_openC(2:end)./msft_closeC(1:end-1))).^2);
msft_prox3C=sqrt(msft_prox3C);
%proxy 4
msft_prox4C=exp(2*log(msft_rangC)-0.86+2*0.29^2);
msft_prox4C=sqrt(msft_prox4C(2:end));


%in sample and forecast proxy
msft_prox1Ci=msft_prox1C(1:nis);msft_prox2Ci=msft_prox2C(1:nis);
msft_prox3Ci=msft_prox3C(1:nis);msft_prox4Ci=msft_prox4C(1:nis);
msft_prox1Cf=msft_prox1C(nis+1:end); msft_prox2Cf=msft_prox2C(nis+1:end); 
msft_prox3Cf=msft_prox3C(nis+1:end); msft_prox4Cf=msft_prox4C(nis+1:end); 

% figure;
% plot(msft_prox1C);hold on;
% plot(msft_prox2C);plot(msft_prox3C);plot(msft_prox3C);

%% ptr proxy
% proxy 1
ptr_prox1C=abs(ptr_rtn-mean(ptr_rtn));
%proxy 2
ptr_highC=PTR(:,3);ptr_lowC=PTR(:,4);ptr_openC=PTR(:,2);ptr_closeC=PTR(:,5);
ptr_rangC=100*log(ptr_highC./ptr_lowC);ptr_rangC(ptr_rangC <= 0)=mean([0 min(ptr_rangC(ptr_rangC>0))]);
ptr_prox2C=sqrt(0.3607*(ptr_rangC.^2));ptr_prox2C=ptr_prox2C(2:end);
%proxy 3
ptr_prox3C=1.107*ptr_prox2C.^2 + 0.68*((100*log(ptr_openC(2:end)./ptr_closeC(1:end-1))).^2);
ptr_prox3C=sqrt(ptr_prox3C);
%proxy 4
ptr_prox4C=exp(2*log(ptr_rangC)-0.86+2*0.29^2);
ptr_prox4C=sqrt(ptr_prox4C(2:end));


%in sample and forecast proxy
ptr_prox1Ci=ptr_prox1C(1:nis);ptr_prox2Ci=ptr_prox2C(1:nis);
ptr_prox3Ci=ptr_prox3C(1:nis);ptr_prox4Ci=ptr_prox4C(1:nis);
ptr_prox1Cf=ptr_prox1C(nis+1:end); ptr_prox2Cf=ptr_prox2C(nis+1:end); 
ptr_prox3Cf=ptr_prox3C(nis+1:end); ptr_prox4Cf=ptr_prox4C(nis+1:end); 
% 
% figure;
% plot(ptr_prox1C);hold on;
% plot(ptr_prox2C);plot(ptr_prox3C);plot(ptr_prox3C);
%% RMSE & MAD Calculation
%% Coke Return
% RMSE
[sqrt(mean((coke_return_forecast_1-coke_rtn_os).^2)) 
sqrt(mean((coke_return_forecast_2-coke_rtn_os).^2)) 
sqrt(mean((coke_return_forecast_4-coke_rtn_os).^2)) 
sqrt(mean((coke_return_forecast_5-coke_rtn_os).^2))]

% MAD
[mean(abs(coke_return_forecast_1-coke_rtn_os)) 
mean(abs(coke_return_forecast_2-coke_rtn_os))
mean(abs(coke_return_forecast_4-coke_rtn_os)) 
mean(abs(coke_return_forecast_5-coke_rtn_os))]

%% Coke Volatility
% Coke Out of Sample forecast RMSE
[sqrt(mean((coke_vol_forecast_1-coke_prox1Cf).^2)) sqrt(mean((coke_vol_forecast_1-coke_prox2Cf).^2)) ...
sqrt(mean((coke_vol_forecast_1-coke_prox3Cf).^2)) sqrt(mean((coke_vol_forecast_1-coke_prox4Cf).^2));
sqrt(mean((coke_vol_forecast_2-coke_prox1Cf).^2)) sqrt(mean((coke_vol_forecast_2-coke_prox2Cf).^2)) ...
sqrt(mean((coke_vol_forecast_2-coke_prox3Cf).^2)) sqrt(mean((coke_vol_forecast_2-coke_prox4Cf).^2));
sqrt(mean((coke_vol_forecast_3-coke_prox1Cf).^2)) sqrt(mean((coke_vol_forecast_3-coke_prox2Cf).^2)) ...
sqrt(mean((coke_vol_forecast_3-coke_prox3Cf).^2)) sqrt(mean((coke_vol_forecast_3-coke_prox4Cf).^2));
sqrt(mean((coke_vol_forecast_4-coke_prox1Cf).^2)) sqrt(mean((coke_vol_forecast_4-coke_prox2Cf).^2)) ...
sqrt(mean((coke_vol_forecast_4-coke_prox3Cf).^2)) sqrt(mean((coke_vol_forecast_4-coke_prox4Cf).^2));
sqrt(mean((coke_vol_forecast_5-coke_prox1Cf).^2)) sqrt(mean((coke_vol_forecast_5-coke_prox2Cf).^2)) ... 
sqrt(mean((coke_vol_forecast_5-coke_prox3Cf).^2)) sqrt(mean((coke_vol_forecast_5-coke_prox4Cf).^2))]

% Coke Out of Sample forecast MAD
[mean(abs(coke_vol_forecast_1-coke_prox1Cf)) mean(abs(coke_vol_forecast_1-coke_prox2Cf)) ...
mean(abs(coke_vol_forecast_1-coke_prox3Cf)) mean(abs(coke_vol_forecast_1-coke_prox4Cf));
mean(abs(coke_vol_forecast_2-coke_prox1Cf)) mean(abs(coke_vol_forecast_2-coke_prox2Cf)) ...
mean(abs(coke_vol_forecast_2-coke_prox3Cf)) mean(abs(coke_vol_forecast_2-coke_prox4Cf));
mean(abs(coke_vol_forecast_3-coke_prox1Cf)) mean(abs(coke_vol_forecast_3-coke_prox2Cf)) ...
mean(abs(coke_vol_forecast_3-coke_prox3Cf)) mean(abs(coke_vol_forecast_3-coke_prox4Cf));
mean(abs(coke_vol_forecast_4-coke_prox1Cf)) mean(abs(coke_vol_forecast_4-coke_prox2Cf)) ...
mean(abs(coke_vol_forecast_4-coke_prox3Cf)) mean(abs(coke_vol_forecast_4-coke_prox4Cf));
mean(abs(coke_vol_forecast_5-coke_prox1Cf)) mean(abs(coke_vol_forecast_5-coke_prox2Cf)) ...
mean(abs(coke_vol_forecast_5-coke_prox3Cf)) mean(abs(coke_vol_forecast_5-coke_prox4Cf))]

%% msft Return
% RMSE
[sqrt(mean((msft_return_forecast_1-msft_rtn_os).^2))
sqrt(mean((msft_return_forecast_2-msft_rtn_os).^2))
sqrt(mean((msft_return_forecast_4-msft_rtn_os).^2)) 
sqrt(mean((msft_return_forecast_5-msft_rtn_os).^2))]

% MAD
[mean(abs(msft_return_forecast_1-msft_rtn_os))
mean(abs(msft_return_forecast_2-msft_rtn_os))
mean(abs(msft_return_forecast_4-msft_rtn_os))
mean(abs(msft_return_forecast_5-msft_rtn_os))]

%% MSFT volatility
% msft Out of Sample forecast RMSE
[sqrt(mean((msft_vol_forecast_1-msft_prox1Cf).^2)) sqrt(mean((msft_vol_forecast_1-msft_prox2Cf).^2)) ...
sqrt(mean((msft_vol_forecast_1-msft_prox3Cf).^2)) sqrt(mean((msft_vol_forecast_1-msft_prox4Cf).^2));
sqrt(mean((msft_vol_forecast_2-msft_prox1Cf).^2)) sqrt(mean((msft_vol_forecast_2-msft_prox2Cf).^2)) ...
sqrt(mean((msft_vol_forecast_2-msft_prox3Cf).^2)) sqrt(mean((msft_vol_forecast_2-msft_prox4Cf).^2));
sqrt(mean((msft_vol_forecast_3-msft_prox1Cf).^2)) sqrt(mean((msft_vol_forecast_3-msft_prox2Cf).^2)) ...
sqrt(mean((msft_vol_forecast_3-msft_prox3Cf).^2)) sqrt(mean((msft_vol_forecast_3-msft_prox4Cf).^2));
sqrt(mean((msft_vol_forecast_4-msft_prox1Cf).^2)) sqrt(mean((msft_vol_forecast_4-msft_prox2Cf).^2)) ...
sqrt(mean((msft_vol_forecast_4-msft_prox3Cf).^2)) sqrt(mean((msft_vol_forecast_4-msft_prox4Cf).^2));
sqrt(mean((msft_vol_forecast_5-msft_prox1Cf).^2)) sqrt(mean((msft_vol_forecast_5-msft_prox2Cf).^2)) ... 
sqrt(mean((msft_vol_forecast_5-msft_prox3Cf).^2)) sqrt(mean((msft_vol_forecast_5-msft_prox4Cf).^2))]

% msft Out of Sample forecast MAD
[mean(abs(msft_vol_forecast_1-msft_prox1Cf)) mean(abs(msft_vol_forecast_1-msft_prox2Cf)) ...
mean(abs(msft_vol_forecast_1-msft_prox3Cf)) mean(abs(msft_vol_forecast_1-msft_prox4Cf));
mean(abs(msft_vol_forecast_2-msft_prox1Cf)) mean(abs(msft_vol_forecast_2-msft_prox2Cf)) ...
mean(abs(msft_vol_forecast_2-msft_prox3Cf)) mean(abs(msft_vol_forecast_2-msft_prox4Cf));
mean(abs(msft_vol_forecast_3-msft_prox1Cf)) mean(abs(msft_vol_forecast_3-msft_prox2Cf)) ...
mean(abs(msft_vol_forecast_3-msft_prox3Cf)) mean(abs(msft_vol_forecast_3-msft_prox4Cf));
mean(abs(msft_vol_forecast_4-msft_prox1Cf)) mean(abs(msft_vol_forecast_4-msft_prox2Cf)) ...
mean(abs(msft_vol_forecast_4-msft_prox3Cf)) mean(abs(msft_vol_forecast_4-msft_prox4Cf));
mean(abs(msft_vol_forecast_5-msft_prox1Cf)) mean(abs(msft_vol_forecast_5-msft_prox2Cf)) ...
mean(abs(msft_vol_forecast_5-msft_prox3Cf)) mean(abs(msft_vol_forecast_5-msft_prox4Cf))]

%% JNJ Return
% RMSE
[sqrt(mean((jnj_return_forecast_1-jnj_rtn_os).^2))
sqrt(mean((jnj_return_forecast_2-jnj_rtn_os).^2))
sqrt(mean((jnj_return_forecast_4-jnj_rtn_os).^2))
sqrt(mean((jnj_return_forecast_5-jnj_rtn_os).^2))]

% MAD
[mean(abs(jnj_return_forecast_1-jnj_rtn_os))
mean(abs(jnj_return_forecast_2-jnj_rtn_os))
mean(abs(jnj_return_forecast_4-jnj_rtn_os))
mean(abs(jnj_return_forecast_5-jnj_rtn_os))]

%% JNJ Volatility
% jnj Out of Sample forecas RMSE
[sqrt(mean((jnj_vol_forecast_1-jnj_prox1Cf).^2)) sqrt(mean((jnj_vol_forecast_1-jnj_prox2Cf).^2)) ...
sqrt(mean((jnj_vol_forecast_1-jnj_prox3Cf).^2)) sqrt(mean((jnj_vol_forecast_1-jnj_prox4Cf).^2));
sqrt(mean((jnj_vol_forecast_2-jnj_prox1Cf).^2)) sqrt(mean((jnj_vol_forecast_2-jnj_prox2Cf).^2)) ...
sqrt(mean((jnj_vol_forecast_2-jnj_prox3Cf).^2)) sqrt(mean((jnj_vol_forecast_2-jnj_prox4Cf).^2));
sqrt(mean((jnj_vol_forecast_3-jnj_prox1Cf).^2)) sqrt(mean((jnj_vol_forecast_3-jnj_prox2Cf).^2)) ...
sqrt(mean((jnj_vol_forecast_3-jnj_prox3Cf).^2)) sqrt(mean((jnj_vol_forecast_3-jnj_prox4Cf).^2));
sqrt(mean((jnj_vol_forecast_4-jnj_prox1Cf).^2)) sqrt(mean((jnj_vol_forecast_4-jnj_prox2Cf).^2)) ...
sqrt(mean((jnj_vol_forecast_4-jnj_prox3Cf).^2)) sqrt(mean((jnj_vol_forecast_4-jnj_prox4Cf).^2));
sqrt(mean((jnj_vol_forecast_5-jnj_prox1Cf).^2)) sqrt(mean((jnj_vol_forecast_5-jnj_prox2Cf).^2)) ... 
sqrt(mean((jnj_vol_forecast_5-jnj_prox3Cf).^2)) sqrt(mean((jnj_vol_forecast_5-jnj_prox4Cf).^2))]

% jnj Out of Sample forecast MAD
[mean(abs(jnj_vol_forecast_1-msjft_prox1Cf)) mean(abs(jnj_vol_forecast_1-jnj_prox2Cf)) ...
mean(abs(jnj_vol_forecast_1-jnj_prox3Cf)) mean(abs(jnj_vol_forecast_1-jnj_prox4Cf));
mean(abs(jnj_vol_forecast_2-jnj_prox1Cf)) mean(abs(jnj_vol_forecast_2-jnj_prox2Cf)) ...
mean(abs(jnj_vol_forecast_2-jnj_prox3Cf)) mean(abs(jnj_vol_forecast_2-jnj_prox4Cf));
mean(abs(jnj_vol_forecast_3-jnj_prox1Cf)) mean(abs(jnj_vol_forecast_3-jnj_prox2Cf)) ...
mean(abs(jnj_vol_forecast_3-jnj_prox3Cf)) mean(abs(jnj_vol_forecast_3-jnj_prox4Cf));
mean(abs(jnj_vol_forecast_4-jnj_prox1Cf)) mean(abs(jnj_vol_forecast_4-jnj_prox2Cf)) ...
mean(abs(jnj_vol_forecast_4-jnj_prox3Cf)) mean(abs(jnj_vol_forecast_4-jnj_prox4Cf));
mean(abs(jnj_vol_forecast_5-jnj_prox1Cf)) mean(abs(jnj_vol_forecast_5-jnj_prox2Cf)) ...
mean(abs(jnj_vol_forecast_5-jnj_prox3Cf)) mean(abs(jnj_vol_forecast_5-jnj_prox4Cf))]

%% mmm Return
% RMSE
[sqrt(mean((mmm_return_forecast_1-mmm_rtn_os).^2))
sqrt(mean((mmm_return_forecast_3-mmm_rtn_os).^2))
sqrt(mean((mmm_return_forecast_4-mmm_rtn_os).^2))
sqrt(mean((mmm_return_forecast_5-mmm_rtn_os).^2))]

% MAD
[mean(abs(mmm_return_forecast_1-mmm_rtn_os))
mean(abs(mmm_return_forecast_3-mmm_rtn_os))
mean(abs(mmm_return_forecast_4-mmm_rtn_os))
mean(abs(mmm_return_forecast_5-mmm_rtn_os))]

%% mmm volatility
[sqrt(mean((mmm_vol_forecast_1-mmm_prox1Cf).^2)) sqrt(mean((mmm_vol_forecast_1-mmm_prox2Cf).^2)) ...
sqrt(mean((mmm_vol_forecast_1-mmm_prox3Cf).^2)) sqrt(mean((mmm_vol_forecast_1-mmm_prox4Cf).^2));
sqrt(mean((mmm_vol_forecast_2-mmm_prox1Cf).^2)) sqrt(mean((mmm_vol_forecast_2-mmm_prox2Cf).^2)) ...
sqrt(mean((mmm_vol_forecast_2-mmm_prox3Cf).^2)) sqrt(mean((mmm_vol_forecast_2-mmm_prox4Cf).^2));
sqrt(mean((mmm_vol_forecast_3-mmm_prox1Cf).^2)) sqrt(mean((mmm_vol_forecast_3-mmm_prox2Cf).^2)) ...
sqrt(mean((mmm_vol_forecast_3-mmm_prox3Cf).^2)) sqrt(mean((mmm_vol_forecast_3-mmm_prox4Cf).^2));
sqrt(mean((mmm_vol_forecast_4-mmm_prox1Cf).^2)) sqrt(mean((mmm_vol_forecast_4-mmm_prox2Cf).^2)) ...
sqrt(mean((mmm_vol_forecast_4-mmm_prox3Cf).^2)) sqrt(mean((mmm_vol_forecast_4-mmm_prox4Cf).^2));
sqrt(mean((mmm_vol_forecast_5-mmm_prox1Cf).^2)) sqrt(mean((mmm_vol_forecast_5-mmm_prox2Cf).^2)) ... 
sqrt(mean((mmm_vol_forecast_5-mmm_prox3Cf).^2)) sqrt(mean((mmm_vol_forecast_5-mmm_prox4Cf).^2))]

% mmm Out of Sample forecast MAD
[mean(abs(mmm_vol_forecast_1-msjft_prox1Cf)) mean(abs(mmm_vol_forecast_1-mmm_prox2Cf)) ...
mean(abs(mmm_vol_forecast_1-mmm_prox3Cf)) mean(abs(mmm_vol_forecast_1-mmm_prox4Cf));
mean(abs(mmm_vol_forecast_2-mmm_prox1Cf)) mean(abs(mmm_vol_forecast_2-mmm_prox2Cf)) ...
mean(abs(mmm_vol_forecast_2-mmm_prox3Cf)) mean(abs(mmm_vol_forecast_2-mmm_prox4Cf));
mean(abs(mmm_vol_forecast_3-mmm_prox1Cf)) mean(abs(mmm_vol_forecast_3-mmm_prox2Cf)) ...
mean(abs(mmm_vol_forecast_3-mmm_prox3Cf)) mean(abs(mmm_vol_forecast_3-mmm_prox4Cf));
mean(abs(mmm_vol_forecast_4-mmm_prox1Cf)) mean(abs(mmm_vol_forecast_4-mmm_prox2Cf)) ...
mean(abs(mmm_vol_forecast_4-mmm_prox3Cf)) mean(abs(mmm_vol_forecast_4-mmm_prox4Cf));
mean(abs(mmm_vol_forecast_5-mmm_prox1Cf)) mean(abs(mmm_vol_forecast_5-mmm_prox2Cf)) ...
mean(abs(mmm_vol_forecast_5-mmm_prox3Cf)) mean(abs(mmm_vol_forecast_5-mmm_prox4Cf))]


%% ptr Return
% RMSE
[sqrt(mean((ptr_return_forecast_1-ptr_rtn_os).^2)) 
sqrt(mean((ptr_return_forecast_2-ptr_rtn_os).^2))
sqrt(mean((ptr_return_forecast_4-ptr_rtn_os).^2)) 
sqrt(mean((ptr_return_forecast_5-ptr_rtn_os).^2))]

% MAD
[mean(abs(ptr_return_forecast_1-ptr_rtn_os))
mean(abs(ptr_return_forecast_2-ptr_rtn_os))
mean(abs(ptr_return_forecast_4-ptr_rtn_os))
mean(abs(ptr_return_forecast_5-ptr_rtn_os)) ]

%% PTR Volatility RMSE
[sqrt(mean((ptr_vol_forecast_1-ptr_prox1Cf).^2)) sqrt(mean((ptr_vol_forecast_1-ptr_prox2Cf).^2)) ...
sqrt(mean((ptr_vol_forecast_1-ptr_prox3Cf).^2)) sqrt(mean((ptr_vol_forecast_1-ptr_prox4Cf).^2));
sqrt(mean((ptr_vol_forecast_2-ptr_prox1Cf).^2)) sqrt(mean((ptr_vol_forecast_2-ptr_prox2Cf).^2)) ...
sqrt(mean((ptr_vol_forecast_2-ptr_prox3Cf).^2)) sqrt(mean((ptr_vol_forecast_2-ptr_prox4Cf).^2));
sqrt(mean((ptr_vol_forecast_3-ptr_prox1Cf).^2)) sqrt(mean((ptr_vol_forecast_3-ptr_prox2Cf).^2)) ...
sqrt(mean((ptr_vol_forecast_3-ptr_prox3Cf).^2)) sqrt(mean((ptr_vol_forecast_3-ptr_prox4Cf).^2));
sqrt(mean((ptr_vol_forecast_4-ptr_prox1Cf).^2)) sqrt(mean((ptr_vol_forecast_4-ptr_prox2Cf).^2)) ...
sqrt(mean((ptr_vol_forecast_4-ptr_prox3Cf).^2)) sqrt(mean((ptr_vol_forecast_4-ptr_prox4Cf).^2));
sqrt(mean((ptr_vol_forecast_5-ptr_prox1Cf).^2)) sqrt(mean((ptr_vol_forecast_5-ptr_prox2Cf).^2)) ... 
sqrt(mean((ptr_vol_forecast_5-ptr_prox3Cf).^2)) sqrt(mean((ptr_vol_forecast_5-ptr_prox4Cf).^2))]

% ptr Out of Sample forecast MAD
[mean(abs(ptr_vol_forecast_1-msjft_prox1Cf)) mean(abs(ptr_vol_forecast_1-ptr_prox2Cf)) ...
mean(abs(ptr_vol_forecast_1-ptr_prox3Cf)) mean(abs(ptr_vol_forecast_1-ptr_prox4Cf));
mean(abs(ptr_vol_forecast_2-ptr_prox1Cf)) mean(abs(ptr_vol_forecast_2-ptr_prox2Cf)) ...
mean(abs(ptr_vol_forecast_2-ptr_prox3Cf)) mean(abs(ptr_vol_forecast_2-ptr_prox4Cf));
mean(abs(ptr_vol_forecast_3-ptr_prox1Cf)) mean(abs(ptr_vol_forecast_3-ptr_prox2Cf)) ...
mean(abs(ptr_vol_forecast_3-ptr_prox3Cf)) mean(abs(ptr_vol_forecast_3-ptr_prox4Cf));
mean(abs(ptr_vol_forecast_4-ptr_prox1Cf)) mean(abs(ptr_vol_forecast_4-ptr_prox2Cf)) ...
mean(abs(ptr_vol_forecast_4-ptr_prox3Cf)) mean(abs(ptr_vol_forecast_4-ptr_prox4Cf));
mean(abs(ptr_vol_forecast_5-ptr_prox1Cf)) mean(abs(ptr_vol_forecast_5-ptr_prox2Cf)) ...
mean(abs(ptr_vol_forecast_5-ptr_prox3Cf)) mean(abs(ptr_vol_forecast_5-ptr_prox4Cf))]
