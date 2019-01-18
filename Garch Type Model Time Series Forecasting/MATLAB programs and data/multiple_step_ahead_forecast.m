%% 1. Accuracy Assessment Proxy Calculation
fs=840;
nis=length(coke_rtn)-fs;
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

%% 2. Forecast & Visualization
%%%%%%%%%%%%%%%%%%%%%%% 1. Coke AR(1)-GARCH(1,1)-t%%%%%%%%%%%%%%%%%%%%
Mdl = arima(1,0,0); Mdl.Variance=garch(1,1); Mdl.Distribution='T';  % refit the oject everytime
[EstMdl,~,~,~] = estimate(Mdl,coke_rtn_is);
[coke_return_forecast_multi_1,~,coke_vol_forecast_multi_1]=forecast(EstMdl,840,'Y0',coke_rtn_is);
[~,coke_vol_is_1,~] = infer(EstMdl,coke_rtn_is);
coke_std_is_1=sqrt(coke_vol_is_1);
coke_std_forecast_multi_1=sqrt(coke_vol_forecast_multi_1);
%%%%%%%%%%%%%%%%%%%% 2. Coke AR(1)-GJR-GARCH(1,1)-t5%%%%%%%%%%%%%%%%%%%
Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T';
[EstMdl,~,~,~] = estimate(Mdl,coke_rtn_is);
[coke_return_forecast_multi_2,~,coke_vol_forecast_multi_2]=forecast(EstMdl,840,'Y0',coke_rtn_is);
[~,coke_vol_is_2,~] = infer(EstMdl,coke_rtn_is);
coke_std_is_2=sqrt(coke_vol_is_2);
coke_std_forecast_multi_2=sqrt(coke_vol_forecast_multi_2);
%%%%%%%%%%%%%%%%%%%% 3. I-GARCH(1,1)-N (For volatility only)%%%%%%%%%%%
[PARAMETERS,~,HT,~,~,~,~]=igarch(coke_rtn_is,1,1,[],[],0);
coke_vol_forecast_multi_3(1)=PARAMETERS*coke_rtn_is(end)^2+(1-PARAMETERS)*HT(end);
coke_vol_forecast_multi_3(2:840)=HT(end);
coke_vol_forecast_multi_3=coke_vol_forecast_multi_3';
coke_std_forecast_multi_3=sqrt(coke_vol_forecast_multi_3);
%%%%%%%%%%%%%%%%%%%%% 4. Ad(25) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coke_return_hs25 = tsmovavg(coke_rtn_is,'w',repelem(1/25, 25),1);
coke_return_forecast_multi_4(1:840,1)=coke_return_hs25(end);
coke_vol_forecast_multi_4(1:840,1)=coke_hs25(end-840);
coke_std_forecast_multi_4=sqrt(coke_vol_forecast_multi_4);
%%%%%%%%%%%%%%%%%%%%%% 5. Ad(5) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coke_return_hs5 = tsmovavg(coke_rtn_is,'w',repelem(1/5,5),1);
coke_return_forecast_multi_5(1:840,1)=coke_return_hs5(end);
coke_vol_forecast_multi_5(1:840,1)=coke_hs5(end-840);
coke_std_forecast_multi_5=sqrt(coke_vol_forecast_multi_5);
% Visualization
% return
figure; plot(date_os(1:10),coke_return_forecast_multi_1(1:10),'r--'); hold on;
plot(date_os(1:10),coke_return_forecast_multi_2(1:10),'k--');
plot(date_os(1:10),coke_return_forecast_multi_4(1:10),'m--');
plot(date_os(1:10),coke_return_forecast_multi_5(1:10),'g--'); 
plot(date_os(1:10),coke_rtn_os(1:10),'y','LineWidth',1.5);
legend('AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T','HS-25','HS-5','Log Return','Location','northeast')
title('Coke Multiple Step Forecast and Log Return');    

% volatility
figure;
plot(date_os,coke_rtn_os,'Color',[0,0.7,0.9]); hold on;
plot(date_os,coke_prox1Cf,'y'); 
plot(date_os,coke_prox2Cf,'Color',[0.9290, 0.6940, 0.1250]); 
plot(date_os,coke_std_forecast_multi_1,'b--','LineWidth',1);
plot(date_os,coke_std_forecast_multi_2,'k--','LineWidth',1);
plot(date_os,coke_std_forecast_multi_3,'c--','LineWidth',1);
plot(date_os,coke_std_forecast_multi_4,'m--','LineWidth',1);
plot(date_os,coke_std_forecast_multi_5,'g--','LineWidth',1); 
legend('Log-Rtn','Proxy1','Proxy2','AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T',...
    'IGARCh(1,1)-G','HS-25','HS-5','Location','southwest')
title('Coke Multiple Step Forecast and std Proxy1&Proxy2');  

% volatility Whole
figure;
plot(date_is(2:end),coke_rtn_is,'Color',[0,0.7,0.9]); hold on;
plot(date_os,coke_rtn_os,'Color',[0,0.7,0.9]);
plot(date_os,coke_prox1Cf,'y');
plot(date_os,coke_prox2Cf,'Color',[0.9290, 0.6940, 0.1250]);
plot(date_os,coke_std_forecast_multi_1,'b--','LineWidth',1);
plot(date_os,coke_std_forecast_multi_2,'k--','LineWidth',1);
plot(date_os,coke_std_forecast_multi_3,'c--','LineWidth',1);
plot(date_os,coke_std_forecast_multi_4,'m--','LineWidth',1);
plot(date_os,coke_std_forecast_multi_5,'g--','LineWidth',1); 
legend('Log-Rtn-IS','Log-Rtn-OS','Proxy2','Proxy1','AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T',...
    'IGARCH(1,1)-G','HS-25','HS-5','Location','northeast')
plot(date_is(2:end),coke_prox1Ci,'y');
plot(date_is(2:end),coke_prox2Ci,'Color',[0.9290, 0.6940, 0.1250]);
plot(date_is(2:end),coke_std_is_1,'b--','LineWidth',1);
plot(date_is(2:end),coke_std_is_2,'k--','LineWidth',1);
plot(date_is(2:end),sqrt(HT),'c--','LineWidth',1);
title('Coke in Sample and Multiple Step std Forecast and std Proxy1&Proxy2'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% MSFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% 1. MSFT AR(1)-GARCH(1,1)-t %%%%%%%%%%%%%%%%%%%%%%
Mdl = arima(1,0,0); Mdl.Variance=garch(1,1); Mdl.Distribution='T';
[EstMdl,~,~,~] = estimate(Mdl,msft_rtn_is);
[msft_return_forecast_multi_1,~,msft_vol_forecast_multi_1]=forecast(EstMdl,840,'Y0',msft_rtn_is);
msft_std_forecast_multi_1=sqrt(msft_vol_forecast_multi_1);
%%%%%%%%%%%%msft_vol_forecast_multi_1%%%%%%%% 2. msft AR(1)-GJR-GARCH(1,1)-t5%%%%%%%%%%%%%%%%%%%
Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T';
[EstMdl,~,~,~] = estimate(Mdl,msft_rtn_is);
[msft_return_forecast_multi_2,~,msft_vol_forecast_multi_2]=forecast(EstMdl,840,'Y0',msft_rtn_is);
msft_std_forecast_multi_2=sqrt(msft_vol_forecast_multi_2);
%%%%%%%%%%%%%%%%%%%% 3. I-GARCH(1,1)-N (For volatility only)%%%%%%%%%%%
[PARAMETERS,~,HT,~,~,~,~]=igarch(msft_rtn_is,1,1,[],[],0);
msft_vol_forecast_multi_3(1)=PARAMETERS*msft_rtn_is(end)^2+(1-PARAMETERS)*HT(end);
msft_vol_forecast_multi_3(2:840)=HT(end);
msft_vol_forecast_multi_3=msft_vol_forecast_multi_3';
msft_std_forecast_multi_3=sqrt(msft_vol_forecast_multi_3);
%%%%%%%%%%%%%%%%%%%%% 4. Ad(25) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msft_return_hs25 = tsmovavg(msft_rtn_is,'w',repelem(1/25, 25),1);
msft_return_forecast_multi_4(1:840,1)=msft_return_hs25(end);
msft_vol_forecast_multi_4(1:840,1)=msft_hs25(end-840);
msft_std_forecast_multi_4=sqrt(msft_vol_forecast_multi_4);
%%%%%%%%%%%%%%%%%%%%%% 5. Ad(5) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msft_return_hs5 = tsmovavg(msft_rtn_is,'w',repelem(1/5,5),1);
msft_return_forecast_multi_5(1:840,1)=msft_return_hs5(end);
msft_vol_forecast_multi_5(1:840,1)=msft_hs5(end-840);
msft_std_forecast_multi_5=sqrt(msft_vol_forecast_multi_5);
% Visualization
% return
figure; plot(date_os(1:10),msft_return_forecast_multi_1(1:10),'r--'); hold on;
plot(date_os(1:10),msft_return_forecast_multi_2(1:10),'k--');
plot(date_os(1:10),msft_return_forecast_multi_4(1:10),'m--');
plot(date_os(1:10),msft_return_forecast_multi_5(1:10),'g--'); 
plot(date_os(1:10),msft_rtn_os(1:10),'y','LineWidth',1.5);
legend('AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T','HS-25','HS-5','Log Return','Location','northeast')
title('MSFT Multiple Step Forecast and Log Return');    

% volatility
figure;plot(date_os,msft_prox2Cf,'y'); hold on;
plot(date_os,msft_vol_forecast_multi_1,'b--','LineWidth',1);
plot(date_os,msft_vol_forecast_multi_2,'k--','LineWidth',1);
plot(date_os,msft_vol_forecast_multi_3,'c--','LineWidth',1);
plot(date_os,msft_vol_forecast_multi_4,'m--','LineWidth',1);
plot(date_os,msft_vol_forecast_multi_5,'g--','LineWidth',1); 
legend('Proxy2','AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T','IGARCh(1,1)-G','HS-25','HS-5','Location','northwest')
title('MSFT Multiple Step Forecast and Vol Proxy2');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%% jnj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% 1. jnj AR(1)-GARCH(1,1)-t %%%%%%%%%%%%%%%%%%%%%%
Mdl = arima(1,0,0); Mdl.Variance=garch(1,1); Mdl.Distribution='T';
[EstMdl,~,~,~] = estimate(Mdl,jnj_rtn_is);
[jnj_return_forecast_multi_1,~,jnj_vol_forecast_multi_1]=forecast(EstMdl,840,'Y0',jnj_rtn_is);
jnj_std_forecast_multi_1=sqrt(jnj_vol_forecast_multi_1);
%%%%%%%%%%%%%%%%%%%% 2. jnj AR(1)-GJR-GARCH(1,1)-t5%%%%%%%%%%%%%%%%%%%
Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T';
[EstMdl,~,~,~] = estimate(Mdl,jnj_rtn_is);
[jnj_return_forecast_multi_2,~,jnj_vol_forecast_multi_2]=forecast(EstMdl,840,'Y0',jnj_rtn_is);
jnj_std_forecast_multi_2=sqrt(jnj_vol_forecast_multi_2);
%%%%%%%%%%%%%%%%%%%% 3. I-GARCH(1,1)-N (For volatility only)%%%%%%%%%%%
[PARAMETERS,~,HT,~,~,~,~]=igarch(jnj_rtn_is,1,1,[],[],0);
jnj_vol_forecast_multi_3(1)=PARAMETERS*jnj_rtn_is(end)^2+(1-PARAMETERS)*HT(end);
jnj_vol_forecast_multi_3(2:840)=HT(end);
jnj_vol_forecast_multi_3=jnj_vol_forecast_multi_3';
jnj_std_forecast_multi_3=sqrt(jnj_vol_forecast_multi_3);
%%%%%%%%%%%%%%%%%%%%% 4. Ad(25) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jnj_return_hs25 = tsmovavg(jnj_rtn_is,'w',repelem(1/25, 25),1);
jnj_return_forecast_multi_4(1:840,1)=jnj_return_hs25(end);
jnj_vol_forecast_multi_4(1:840,1)=jnj_hs25(end-840);
jnj_std_forecast_multi_4=sqrt(jnj_vol_forecast_multi_4);
%%%%%%%%%%%%%%%%%%%%%% 5. Ad(5) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jnj_return_hs5 = tsmovavg(jnj_rtn_is,'w',repelem(1/5,5),1);
jnj_return_forecast_multi_5(1:840,1)=jnj_return_hs5(end);
jnj_vol_forecast_multi_5(1:840,1)=jnj_hs5(end-840);
jnj_std_forecast_multi_5=sqrt(jnj_vol_forecast_multi_5);
% Visualization
% return
figure; plot(date_os(1:10),jnj_return_forecast_multi_1(1:10),'r--'); hold on;
plot(date_os(1:10),jnj_return_forecast_multi_2(1:10),'k--');
plot(date_os(1:10),jnj_return_forecast_multi_4(1:10),'m--');
plot(date_os(1:10),jnj_return_forecast_multi_5(1:10),'g--'); 
plot(date_os(1:10),jnj_rtn_os(1:10),'y','LineWidth',1.5);
legend('AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T','HS-25','HS-5','Log Return','Location','northeast')
title('JNJ Multiple Step Forecast and Log Return');    

% volatility
figure;plot(date_os,jnj_prox2Cf,'y'); hold on;
plot(date_os,jnj_vol_forecast_multi_1,'b--','LineWidth',1);
plot(date_os,jnj_vol_forecast_multi_2,'k--','LineWidth',1);
plot(date_os,jnj_vol_forecast_multi_3,'c--','LineWidth',1);
plot(date_os,jnj_vol_forecast_multi_4,'m--','LineWidth',1);
plot(date_os,jnj_vol_forecast_multi_5,'g--','LineWidth',1); 
legend('Proxy2','AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T','IGARCh(1,1)-G','HS-25','HS-5','Location','northwest')
title('JNJ Multiple Step Forecast and Vol Proxy2');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%% mmm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% 1. mmm AR(1)-GARCH(1,1)-t %%%%%%%%%%%%%%%%%%%%%%
Mdl = arima(1,0,0); Mdl.Variance=garch(1,1); Mdl.Distribution='T';
[EstMdl,~,~,~] = estimate(Mdl,mmm_rtn_is);
[mmm_return_forecast_multi_1,~,mmm_vol_forecast_multi_1]=forecast(EstMdl,840,'Y0',mmm_rtn_is);
mmm_std_forecast_multi_1=sqrt(mmm_vol_forecast_multi_1);
% Visualization
% return
figure; plot(date_os(1:10),mmm_return_forecast_multi_1(1:10),'r--'); hold on;
plot(date_os(1:10),mmm_return_forecast_multi_2(1:10),'k--');
plot(date_os(1:10),mmm_return_forecast_multi_4(1:10),'m--');
plot(date_os(1:10),mmm_return_forecast_multi_5(1:10),'g--'); 
plot(date_os(1:10),mmm_rtn_os(1:10),'y','LineWidth',1.5);
legend('AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T','HS-25','HS-5','Log Return','Location','northeast')
title('mmm Multiple Step Forecast and Log Return');    

% volatility
figure;plot(date_os,mmm_prox2Cf,'y'); hold on;
plot(date_os,mmm_vol_forecast_multi_1,'b--','LineWidth',1);
plot(date_os,mmm_vol_forecast_multi_2,'k--','LineWidth',1);
plot(date_os,mmm_vol_forecast_multi_3,'c--','LineWidth',1);
plot(date_os,mmm_vol_forecast_multi_4,'m--','LineWidth',1);
plot(date_os,mmm_vol_forecast_multi_5,'g--','LineWidth',1); 
legend('Proxy2','AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T','IGARCh(1,1)-G','HS-25','HS-5','Location','northwest')
title('mmm Multiple Step Forecast and Vol Proxy2');  

%%%%%%%%%%%%%%%%%%%% 2. mmm AR(1)-GJR-GARCH(1,1)-t5%%%%%%%%%%%%%%%%%%%
Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T';
[EstMdl,~,~,~] = estimate(Mdl,mmm_rtn_is);
[mmm_return_forecast_multi_2,~,mmm_vol_forecast_multi_2]=forecast(EstMdl,840,'Y0',mmm_rtn_is);
mmm_std_forecast_multi_2=sqrt(mmm_vol_forecast_multi_2);
%%%%%%%%%%%%%%%%%%%% 3. I-GARCH(1,1)-N (For volatility only)%%%%%%%%%%%
[PARAMETERS,~,HT,~,~,~,~]=igarch(mmm_rtn_is,1,1,[],[],0);
mmm_vol_forecast_multi_3(1)=PARAMETERS*mmm_rtn_is(end)^2+(1-PARAMETERS)*HT(end);
mmm_vol_forecast_multi_3(2:840)=HT(end);
mmm_vol_forecast_multi_3=mmm_vol_forecast_multi_3';
mmm_std_forecast_multi_3=sqrt(mmm_vol_forecast_multi_3);
%%%%%%%%%%%%%%%%%%%%% 4. Ad(25) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmm_return_hs25 = tsmovavg(mmm_rtn_is,'w',repelem(1/25, 25),1);
mmm_return_forecast_multi_4(1:840,1)=mmm_return_hs25(end);
mmm_vol_forecast_multi_4(1:840,1)=mmm_hs25(end-840);
mmm_std_forecast_multi_4=sqrt(mmm_vol_forecast_multi_4);
%%%%%%%%%%%%%%%%%%%%%% 5. Ad(5) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mmm_return_hs5 = tsmovavg(mmm_rtn_is,'w',repelem(1/5,5),1);
mmm_return_forecast_multi_5(1:840,1)=mmm_return_hs5(end);
mmm_vol_forecast_multi_5(1:840,1)=mmm_hs5(end-840);    
mmm_std_forecast_multi_5=sqrt(mmm_vol_forecast_multi_5);
% Visualization
% return
figure; plot(date_os(1:10),mmm_return_forecast_multi_1(1:10),'r--'); hold on;
plot(date_os(1:10),mmm_return_forecast_multi_2(1:10),'k--');
plot(date_os(1:10),mmm_return_forecast_multi_4(1:10),'m--');
plot(date_os(1:10),mmm_return_forecast_multi_5(1:10),'g--'); 
plot(date_os(1:10),mmm_rtn_os(1:10),'y','LineWidth',1.5);
legend('AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T','HS-25','HS-5','Log Return','Location','northeast')
title('MMM Multiple Step Forecast and Log Return');    

% volatility
figure;plot(date_os,mmm_prox2Cf,'y'); hold on;
plot(date_os,mmm_vol_forecast_multi_1,'b--','LineWidth',1);
plot(date_os,mmm_vol_forecast_multi_2,'k--','LineWidth',1);
plot(date_os,mmm_vol_forecast_multi_3,'c--','LineWidth',1);
plot(date_os,mmm_vol_forecast_multi_4,'m--','LineWidth',1);
plot(date_os,mmm_vol_forecast_multi_5,'g--','LineWidth',1); 
legend('Proxy2','AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T','IGARCh(1,1)-G','HS-25','HS-5','Location','northwest')
title('MMM Multiple Step Forecast and Vol Proxy2');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%% ptr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% 1. ptr AR(1)-GARCH(1,1)-t %%%%%%%%%%%%%%%%%%%%%%
Mdl = arima(1,0,0); Mdl.Variance=garch(1,1); Mdl.Distribution='T';
[EstMdl,~,~,~] = estimate(Mdl,ptr_rtn_is);
[ptr_return_forecast_multi_1,~,ptr_vol_forecast_multi_1]=forecast(EstMdl,840,'Y0',ptr_rtn_is);
[~,ptr_vol_is_1,~] = infer(EstMdl,ptr_rtn_is);
ptr_std_forecast_multi_1=sqrt(ptr_vol_forecast_multi_1);
%%%%%%%%%%%%%%%%%%%% 2. ptr AR(1)-GJR-GARCH(1,1)-t5%%%%%%%%%%%%%%%%%%%
Mdl = arima(1,0,0); Mdl.Variance=gjr(1,1); Mdl.Distribution='T';
[EstMdl,~,~,~] = estimate(Mdl,ptr_rtn_is);
[ptr_return_forecast_multi_2,~,ptr_vol_forecast_multi_2]=forecast(EstMdl,840,'Y0',ptr_rtn_is);
[~,ptr_vol_is_2,~] = infer(EstMdl,ptr_rtn_is);
ptr_std_forecast_multi_2=sqrt(ptr_vol_forecast_multi_2);
%%%%%%%%%%%%%%%%%%%% 3. I-GARCH(1,1)-N (For volatility only)%%%%%%%%%%%
[PARAMETERS,~,HT,~,~,~,~]=igarch(ptr_rtn_is,1,1,[],[],0);
ptr_vol_forecast_multi_3(1)=PARAMETERS*ptr_rtn_is(end)^2+(1-PARAMETERS)*HT(end);
ptr_vol_forecast_multi_3(2:840)=HT(end);
ptr_vol_forecast_multi_3=ptr_vol_forecast_multi_3';
ptr_std_forecast_multi_3=sqrt(ptr_vol_forecast_multi_3);
%%%%%%%%%%%%%%%%%%%%% 4. Ad(25) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ptr_return_hs25 = tsmovavg(ptr_rtn_is,'w',repelem(1/25, 25),1);
ptr_return_forecast_multi_4(1:840,1)=ptr_return_hs25(end);
ptr_vol_forecast_multi_4(1:840,1)=ptr_hs25(end-840);
ptr_std_forecast_multi_4=sqrt(ptr_vol_forecast_multi_4);
%%%%%%%%%%%%%%%%%%%%%% 5. Ad(5) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ptr_return_hs5 = tsmovavg(ptr_rtn_is,'w',repelem(1/5,5),1);
ptr_return_forecast_multi_5(1:840,1)=ptr_return_hs5(end);
ptr_vol_forecast_multi_5(1:840,1)=ptr_hs5(end-840); 
ptr_std_forecast_multi_5=sqrt(ptr_vol_forecast_multi_5);
% Visualization
% return
figure; plot(date_os(1:10),ptr_return_forecast_multi_1(1:10),'r--'); hold on;
plot(date_os(1:10),ptr_return_forecast_multi_2(1:10),'k--');
plot(date_os(1:10),ptr_return_forecast_multi_4(1:10),'m--');
plot(date_os(1:10),ptr_return_forecast_multi_5(1:10),'g--'); 
plot(date_os(1:10),ptr_rtn_os(1:10),'y','LineWidth',1.5);
legend('AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T','HS-25','HS-5','Log Return','Location','northeast')
title('ptr Multiple Step Forecast and Log Return');    

% volatility
plot(date_os,ptr_prox1Cf,'y');  hold on;
plot(date_os,ptr_prox2Cf,'Color',[0.9290, 0.6940, 0.1250]);
plot(date_os,ptr_vol_forecast_multi_1,'b--','LineWidth',1);
plot(date_os,ptr_vol_forecast_multi_2,'k--','LineWidth',1);
plot(date_os,ptr_vol_forecast_multi_3,'c--','LineWidth',1);
plot(date_os,ptr_vol_forecast_multi_4,'m--','LineWidth',1);
plot(date_os,ptr_vol_forecast_multi_5,'g--','LineWidth',1); 
legend('Proxy2','Proxy1','AR(1)-GARCH(1,1)-T','AR(1)-GJR-GARCH(1,1)-T',...
    'IGARCH(1,1)-G','HS-25','HS-5','Location','northwest')
plot(date_is(2:end),ptr_prox1Ci,'y');
plot(date_is(2:end),ptr_prox2Ci,'Color',[0.9290, 0.6940, 0.1250]);
plot(date_is(2:end),ptr_vol_is_1,'b--','LineWidth',1);
plot(date_is(2:end),ptr_vol_is_2,'k--','LineWidth',1);
plot(date_is(2:end),HT,'c--','LineWidth',1);
title('ptr In Sample and Multiple Step Vol Forecast and Vol Proxy1&Proxy2');  
%% Accuracy Assessment
%% coke
% Rtn
% rmse_multi
rtn_rmse_multi(:,1)=[sqrt(mean((coke_return_forecast_multi_1-coke_rtn_os).^2)) 
sqrt(mean((coke_return_forecast_multi_2-coke_rtn_os).^2)) 
sqrt(mean((coke_return_forecast_multi_4-coke_rtn_os).^2)) 
sqrt(mean((coke_return_forecast_multi_5-coke_rtn_os).^2))]

% mad_multi
rtn_mad_multi(:,1)=[mean(abs(coke_return_forecast_multi_1-coke_rtn_os)) 
mean(abs(coke_return_forecast_multi_2-coke_rtn_os))
mean(abs(coke_return_forecast_multi_4-coke_rtn_os)) 
mean(abs(coke_return_forecast_multi_5-coke_rtn_os))]
% Vol
% Proxy1
std_rmse_multi_prox1(:,1)=[sqrt(mean((coke_std_forecast_multi_1-coke_prox1Cf).^2)) 
              sqrt(mean((coke_std_forecast_multi_2-coke_prox1Cf).^2))
              sqrt(mean((coke_std_forecast_multi_3-coke_prox1Cf).^2))
              sqrt(mean((coke_std_forecast_multi_4-coke_prox1Cf).^2))
              sqrt(mean((coke_std_forecast_multi_5-coke_prox1Cf).^2))]
% Proxy2
std_rmse_multi_prox2(:,1)=[sqrt(mean((coke_std_forecast_multi_1-coke_prox2Cf).^2)) 
              sqrt(mean((coke_std_forecast_multi_2-coke_prox2Cf).^2))
              sqrt(mean((coke_std_forecast_multi_3-coke_prox2Cf).^2))
              sqrt(mean((coke_std_forecast_multi_4-coke_prox2Cf).^2))
              sqrt(mean((coke_std_forecast_multi_5-coke_prox2Cf).^2))]
% Proxy3
std_rmse_multi_prox3(:,1)=[sqrt(mean((coke_std_forecast_multi_1-coke_prox3Cf).^2)) 
              sqrt(mean((coke_std_forecast_multi_2-coke_prox3Cf).^2))
              sqrt(mean((coke_std_forecast_multi_3-coke_prox3Cf).^2))
              sqrt(mean((coke_std_forecast_multi_4-coke_prox3Cf).^2))
              sqrt(mean((coke_std_forecast_multi_5-coke_prox3Cf).^2))]
% Proxy4
std_rmse_multi_prox4(:,1)=[sqrt(mean((coke_std_forecast_multi_1-coke_prox4Cf).^2)) 
              sqrt(mean((coke_std_forecast_multi_2-coke_prox4Cf).^2))
              sqrt(mean((coke_std_forecast_multi_3-coke_prox4Cf).^2))
              sqrt(mean((coke_std_forecast_multi_4-coke_prox4Cf).^2))
              sqrt(mean((coke_std_forecast_multi_5-coke_prox4Cf).^2))]
          
% mad_multi
% Proxy1
std_mad_multi_prox1(:,1)=[mean(abs(coke_std_forecast_multi_1-coke_prox1Cf))
                    mean(abs(coke_std_forecast_multi_2-coke_prox1Cf))
                    mean(abs(coke_std_forecast_multi_3-coke_prox1Cf))
                    mean(abs(coke_std_forecast_multi_4-coke_prox1Cf))
                    mean(abs(coke_std_forecast_multi_5-coke_prox1Cf))]
% Proxy2
std_mad_multi_prox2(:,1)=[mean(abs(coke_std_forecast_multi_1-coke_prox2Cf))
                    mean(abs(coke_std_forecast_multi_2-coke_prox2Cf))
                    mean(abs(coke_std_forecast_multi_3-coke_prox2Cf))
                    mean(abs(coke_std_forecast_multi_4-coke_prox2Cf))
                    mean(abs(coke_std_forecast_multi_5-coke_prox2Cf))]
% Proxy3
std_mad_multi_prox3(:,1)=[mean(abs(coke_std_forecast_multi_1-coke_prox3Cf))
                    mean(abs(coke_std_forecast_multi_2-coke_prox3Cf))
                    mean(abs(coke_std_forecast_multi_3-coke_prox3Cf))
                    mean(abs(coke_std_forecast_multi_4-coke_prox3Cf))
                    mean(abs(coke_std_forecast_multi_5-coke_prox3Cf))]
% Proxy4
std_mad_multi_prox4(:,1)=[mean(abs(coke_std_forecast_multi_1-coke_prox4Cf))
                    mean(abs(coke_std_forecast_multi_2-coke_prox4Cf))
                    mean(abs(coke_std_forecast_multi_3-coke_prox4Cf))
                    mean(abs(coke_std_forecast_multi_4-coke_prox4Cf))
                    mean(abs(coke_std_forecast_multi_5-coke_prox4Cf))]   

%% MSFT
% RTN
% rmse_multi
rtn_rmse_multi(:,2)=[sqrt(mean((msft_return_forecast_multi_1-msft_rtn_os).^2))
sqrt(mean((msft_return_forecast_multi_2-msft_rtn_os).^2))
sqrt(mean((msft_return_forecast_multi_4-msft_rtn_os).^2)) 
sqrt(mean((msft_return_forecast_multi_5-msft_rtn_os).^2))]

% mad_multi
rtn_mad_multi(:,2)=[mean(abs(msft_return_forecast_multi_1-msft_rtn_os))
mean(abs(msft_return_forecast_multi_2-msft_rtn_os))
mean(abs(msft_return_forecast_multi_4-msft_rtn_os))
mean(abs(msft_return_forecast_multi_5-msft_rtn_os))]

% VOL
% vol
% rmse_multi
% Proxy1
std_rmse_multi_prox1(:,2)=[sqrt(mean((msft_std_forecast_multi_1-msft_prox1Cf).^2)) 
              sqrt(mean((msft_std_forecast_multi_2-msft_prox1Cf).^2))
              sqrt(mean((msft_std_forecast_multi_3-msft_prox1Cf).^2))
              sqrt(mean((msft_std_forecast_multi_4-msft_prox1Cf).^2))
              sqrt(mean((msft_std_forecast_multi_5-msft_prox1Cf).^2))]
% Proxy2
std_rmse_multi_prox2(:,2)=[sqrt(mean((msft_std_forecast_multi_1-msft_prox2Cf).^2)) 
              sqrt(mean((msft_std_forecast_multi_2-msft_prox2Cf).^2))
              sqrt(mean((msft_std_forecast_multi_3-msft_prox2Cf).^2))
              sqrt(mean((msft_std_forecast_multi_4-msft_prox2Cf).^2))
              sqrt(mean((msft_std_forecast_multi_5-msft_prox2Cf).^2))]
% Proxy3
std_rmse_multi_prox3(:,2)=[sqrt(mean((msft_std_forecast_multi_1-msft_prox3Cf).^2)) 
              sqrt(mean((msft_std_forecast_multi_2-msft_prox3Cf).^2))
              sqrt(mean((msft_std_forecast_multi_3-msft_prox3Cf).^2))
              sqrt(mean((msft_std_forecast_multi_4-msft_prox3Cf).^2))
              sqrt(mean((msft_std_forecast_multi_5-msft_prox3Cf).^2))]
% Proxy4
std_rmse_multi_prox4(:,2)=[sqrt(mean((msft_std_forecast_multi_1-msft_prox4Cf).^2)) 
              sqrt(mean((msft_std_forecast_multi_2-msft_prox4Cf).^2))
              sqrt(mean((msft_std_forecast_multi_3-msft_prox4Cf).^2))
              sqrt(mean((msft_std_forecast_multi_4-msft_prox4Cf).^2))
              sqrt(mean((msft_std_forecast_multi_5-msft_prox4Cf).^2))]

% mad_multi
% Proxy1
std_mad_multi_prox1(:,2)=[mean(abs(msft_std_forecast_multi_1-msft_prox1Cf))
                    mean(abs(msft_std_forecast_multi_2-msft_prox1Cf))
                    mean(abs(msft_std_forecast_multi_3-msft_prox1Cf))
                    mean(abs(msft_std_forecast_multi_4-msft_prox1Cf))
                    mean(abs(msft_std_forecast_multi_5-msft_prox1Cf))]
% Proxy2
std_mad_multi_prox2(:,2)=[mean(abs(msft_std_forecast_multi_1-msft_prox2Cf))
                    mean(abs(msft_std_forecast_multi_2-msft_prox2Cf))
                    mean(abs(msft_std_forecast_multi_3-msft_prox2Cf))
                    mean(abs(msft_std_forecast_multi_4-msft_prox2Cf))
                    mean(abs(msft_std_forecast_multi_5-msft_prox2Cf))]
% Proxy3
std_mad_multi_prox3(:,2)=[mean(abs(msft_std_forecast_multi_1-msft_prox3Cf))
                    mean(abs(msft_std_forecast_multi_2-msft_prox3Cf))
                    mean(abs(msft_std_forecast_multi_3-msft_prox3Cf))
                    mean(abs(msft_std_forecast_multi_4-msft_prox3Cf))
                    mean(abs(msft_std_forecast_multi_5-msft_prox3Cf))]
% Proxy4
std_mad_multi_prox4(:,2)=[mean(abs(msft_std_forecast_multi_1-msft_prox4Cf))
                    mean(abs(msft_std_forecast_multi_2-msft_prox4Cf))
                    mean(abs(msft_std_forecast_multi_3-msft_prox4Cf))
                    mean(abs(msft_std_forecast_multi_4-msft_prox4Cf))
                    mean(abs(msft_std_forecast_multi_5-msft_prox4Cf))]

%% JNJ
% RTN
% rmse_multi
rtn_rmse_multi(:,3)=[sqrt(mean((jnj_return_forecast_multi_1-jnj_rtn_os).^2))
sqrt(mean((jnj_return_forecast_multi_2-jnj_rtn_os).^2))
sqrt(mean((jnj_return_forecast_multi_4-jnj_rtn_os).^2))
sqrt(mean((jnj_return_forecast_multi_5-jnj_rtn_os).^2))]

% mad_multi
rtn_mad_multi(:,3)=[mean(abs(jnj_return_forecast_multi_1-jnj_rtn_os))
mean(abs(jnj_return_forecast_multi_2-jnj_rtn_os))
mean(abs(jnj_return_forecast_multi_4-jnj_rtn_os))
mean(abs(jnj_return_forecast_multi_5-jnj_rtn_os))]

% rmse_multi
% Proxy1
std_rmse_multi_prox1(:,3)=[sqrt(mean((jnj_std_forecast_multi_1-jnj_prox1Cf).^2)) 
              sqrt(mean((jnj_std_forecast_multi_2-jnj_prox1Cf).^2))
              sqrt(mean((jnj_std_forecast_multi_3-jnj_prox1Cf).^2))
              sqrt(mean((jnj_std_forecast_multi_4-jnj_prox1Cf).^2))
              sqrt(mean((jnj_std_forecast_multi_5-jnj_prox1Cf).^2))]
% Proxy2
std_rmse_multi_prox2(:,3)=[sqrt(mean((jnj_std_forecast_multi_1-jnj_prox2Cf).^2)) 
              sqrt(mean((jnj_std_forecast_multi_2-jnj_prox2Cf).^2))
              sqrt(mean((jnj_std_forecast_multi_3-jnj_prox2Cf).^2))
              sqrt(mean((jnj_std_forecast_multi_4-jnj_prox2Cf).^2))
              sqrt(mean((jnj_std_forecast_multi_5-jnj_prox2Cf).^2))]
% Proxy3
std_rmse_multi_prox3(:,3)=[sqrt(mean((jnj_std_forecast_multi_1-jnj_prox3Cf).^2)) 
              sqrt(mean((jnj_std_forecast_multi_2-jnj_prox3Cf).^2))
              sqrt(mean((jnj_std_forecast_multi_3-jnj_prox3Cf).^2))
              sqrt(mean((jnj_std_forecast_multi_4-jnj_prox3Cf).^2))
              sqrt(mean((jnj_std_forecast_multi_5-jnj_prox3Cf).^2))]
% Proxy4
std_rmse_multi_prox4(:,3)=[sqrt(mean((jnj_std_forecast_multi_1-jnj_prox4Cf).^2)) 
              sqrt(mean((jnj_std_forecast_multi_2-jnj_prox4Cf).^2))
              sqrt(mean((jnj_std_forecast_multi_3-jnj_prox4Cf).^2))
              sqrt(mean((jnj_std_forecast_multi_4-jnj_prox4Cf).^2))
              sqrt(mean((jnj_std_forecast_multi_5-jnj_prox4Cf).^2))]
          
% mad_multi
% Proxy1
std_mad_multi_prox1(:,3)=[mean(abs(jnj_std_forecast_multi_1-jnj_prox1Cf))
                    mean(abs(jnj_std_forecast_multi_2-jnj_prox1Cf))
                    mean(abs(jnj_std_forecast_multi_3-jnj_prox1Cf))
                    mean(abs(jnj_std_forecast_multi_4-jnj_prox1Cf))
                    mean(abs(jnj_std_forecast_multi_5-jnj_prox1Cf))]
% Proxy2
std_mad_multi_prox2(:,3)=[mean(abs(jnj_std_forecast_multi_1-jnj_prox2Cf))
                    mean(abs(jnj_std_forecast_multi_2-jnj_prox2Cf))
                    mean(abs(jnj_std_forecast_multi_3-jnj_prox2Cf))
                    mean(abs(jnj_std_forecast_multi_4-jnj_prox2Cf))
                    mean(abs(jnj_std_forecast_multi_5-jnj_prox2Cf))]
% Proxy3
std_mad_multi_prox3(:,3)=[mean(abs(jnj_std_forecast_multi_1-jnj_prox3Cf))
                    mean(abs(jnj_std_forecast_multi_2-jnj_prox3Cf))
                    mean(abs(jnj_std_forecast_multi_3-jnj_prox3Cf))
                    mean(abs(jnj_std_forecast_multi_4-jnj_prox3Cf))
                    mean(abs(jnj_std_forecast_multi_5-jnj_prox3Cf))]
% Proxy4
std_mad_multi_prox4(:,3)=[mean(abs(jnj_std_forecast_multi_1-jnj_prox4Cf))
                    mean(abs(jnj_std_forecast_multi_2-jnj_prox4Cf))
                    mean(abs(jnj_std_forecast_multi_3-jnj_prox4Cf))
                    mean(abs(jnj_std_forecast_multi_4-jnj_prox4Cf))
                    mean(abs(jnj_std_forecast_multi_5-jnj_prox4Cf))]
%% MMM
% RTN
%rmse_multi
rtn_rmse_multi(:,4)=[sqrt(mean((mmm_return_forecast_multi_1-mmm_rtn_os).^2))
sqrt(mean((mmm_return_forecast_multi_2-mmm_rtn_os).^2))
sqrt(mean((mmm_return_forecast_multi_4-mmm_rtn_os).^2))
sqrt(mean((mmm_return_forecast_multi_5-mmm_rtn_os).^2))]

% mad_multi
rtn_mad_multi(:,4)=[mean(abs(mmm_return_forecast_multi_1-mmm_rtn_os))
mean(abs(mmm_return_forecast_multi_2-mmm_rtn_os))
mean(abs(mmm_return_forecast_multi_4-mmm_rtn_os))
mean(abs(mmm_return_forecast_multi_5-mmm_rtn_os))]

% VOL
% rmse_multi
% Proxy1
std_rmse_multi_prox1(:,4)=[sqrt(mean((mmm_std_forecast_multi_1-mmm_prox1Cf).^2)) 
              sqrt(mean((mmm_std_forecast_multi_2-mmm_prox1Cf).^2))
              sqrt(mean((mmm_std_forecast_multi_3-mmm_prox1Cf).^2))
              sqrt(mean((mmm_std_forecast_multi_4-mmm_prox1Cf).^2))
              sqrt(mean((mmm_std_forecast_multi_5-mmm_prox1Cf).^2))]
% Proxy2
std_rmse_multi_prox2(:,4)=[sqrt(mean((mmm_std_forecast_multi_1-mmm_prox2Cf).^2)) 
              sqrt(mean((mmm_std_forecast_multi_2-mmm_prox2Cf).^2))
              sqrt(mean((mmm_std_forecast_multi_3-mmm_prox2Cf).^2))
              sqrt(mean((mmm_std_forecast_multi_4-mmm_prox2Cf).^2))
              sqrt(mean((mmm_std_forecast_multi_5-mmm_prox2Cf).^2))]
% Proxy3
std_rmse_multi_prox3(:,4)=[sqrt(mean((mmm_std_forecast_multi_1-mmm_prox3Cf).^2)) 
              sqrt(mean((mmm_std_forecast_multi_2-mmm_prox3Cf).^2))
              sqrt(mean((mmm_std_forecast_multi_3-mmm_prox3Cf).^2))
              sqrt(mean((mmm_std_forecast_multi_4-mmm_prox3Cf).^2))
              sqrt(mean((mmm_std_forecast_multi_5-mmm_prox3Cf).^2))]
% Proxy4
std_rmse_multi_prox4(:,4)=[sqrt(mean((mmm_std_forecast_multi_1-mmm_prox4Cf).^2)) 
              sqrt(mean((mmm_std_forecast_multi_2-mmm_prox4Cf).^2))
              sqrt(mean((mmm_std_forecast_multi_3-mmm_prox4Cf).^2))
              sqrt(mean((mmm_std_forecast_multi_4-mmm_prox4Cf).^2))
              sqrt(mean((mmm_std_forecast_multi_5-mmm_prox4Cf).^2))]
          
% mad_multi
% Proxy1
std_mad_multi_prox1(:,4)=[mean(abs(mmm_std_forecast_multi_1-mmm_prox1Cf))
                    mean(abs(mmm_std_forecast_multi_2-mmm_prox1Cf))
                    mean(abs(mmm_std_forecast_multi_3-mmm_prox1Cf))
                    mean(abs(mmm_std_forecast_multi_4-mmm_prox1Cf))
                    mean(abs(mmm_std_forecast_multi_5-mmm_prox1Cf))]
% Proxy2
std_mad_multi_prox2(:,4)=[mean(abs(mmm_std_forecast_multi_1-mmm_prox2Cf))
                    mean(abs(mmm_std_forecast_multi_2-mmm_prox2Cf))
                    mean(abs(mmm_std_forecast_multi_3-mmm_prox2Cf))
                    mean(abs(mmm_std_forecast_multi_4-mmm_prox2Cf))
                    mean(abs(mmm_std_forecast_multi_5-mmm_prox2Cf))]
% Proxy3
std_mad_multi_prox3(:,4)=[mean(abs(mmm_std_forecast_multi_1-mmm_prox3Cf))
                    mean(abs(mmm_std_forecast_multi_2-mmm_prox3Cf))
                    mean(abs(mmm_std_forecast_multi_3-mmm_prox3Cf))
                    mean(abs(mmm_std_forecast_multi_4-mmm_prox3Cf))
                    mean(abs(mmm_std_forecast_multi_5-mmm_prox3Cf))]
% Proxy4
std_mad_multi_prox4(:,4)=[mean(abs(mmm_std_forecast_multi_1-mmm_prox4Cf))
                    mean(abs(mmm_std_forecast_multi_2-mmm_prox4Cf))
                    mean(abs(mmm_std_forecast_multi_3-mmm_prox4Cf))
                    mean(abs(mmm_std_forecast_multi_4-mmm_prox4Cf))
                    mean(abs(mmm_std_forecast_multi_5-mmm_prox4Cf))]
                
%% PTR
% RTN
rtn_rmse_multi(:,5)=[sqrt(mean((ptr_return_forecast_multi_1-ptr_rtn_os).^2)) 
sqrt(mean((ptr_return_forecast_multi_2-ptr_rtn_os).^2))
sqrt(mean((ptr_return_forecast_multi_4-ptr_rtn_os).^2)) 
sqrt(mean((ptr_return_forecast_multi_5-ptr_rtn_os).^2))]

% mad_multi
rtn_mad_multi(:,5)=[mean(abs(ptr_return_forecast_multi_1-ptr_rtn_os))
mean(abs(ptr_return_forecast_multi_2-ptr_rtn_os))
mean(abs(ptr_return_forecast_multi_4-ptr_rtn_os))
mean(abs(ptr_return_forecast_multi_5-ptr_rtn_os))]

% VOL
% rmse_multi
% Proxy1
std_rmse_multi_prox1(:,5)=[sqrt(mean((ptr_std_forecast_multi_1-ptr_prox1Cf).^2)) 
              sqrt(mean((ptr_std_forecast_multi_2-ptr_prox1Cf).^2))
              sqrt(mean((ptr_std_forecast_multi_3-ptr_prox1Cf).^2))
              sqrt(mean((ptr_std_forecast_multi_4-ptr_prox1Cf).^2))
              sqrt(mean((ptr_std_forecast_multi_5-ptr_prox1Cf).^2))]
% Proxy2
std_rmse_multi_prox2(:,5)=[sqrt(mean((ptr_std_forecast_multi_1-ptr_prox2Cf).^2)) 
              sqrt(mean((ptr_std_forecast_multi_2-ptr_prox2Cf).^2))
              sqrt(mean((ptr_std_forecast_multi_3-ptr_prox2Cf).^2))
              sqrt(mean((ptr_std_forecast_multi_4-ptr_prox2Cf).^2))
              sqrt(mean((ptr_std_forecast_multi_5-ptr_prox2Cf).^2))]
% Proxy3
std_rmse_multi_prox3(:,5)=[sqrt(mean((ptr_std_forecast_multi_1-ptr_prox3Cf).^2)) 
              sqrt(mean((ptr_std_forecast_multi_2-ptr_prox3Cf).^2))
              sqrt(mean((ptr_std_forecast_multi_3-ptr_prox3Cf).^2))
              sqrt(mean((ptr_std_forecast_multi_4-ptr_prox3Cf).^2))
              sqrt(mean((ptr_std_forecast_multi_5-ptr_prox3Cf).^2))]
% Proxy4
std_rmse_multi_prox4(:,5)=[sqrt(mean((ptr_std_forecast_multi_1-ptr_prox4Cf).^2)) 
              sqrt(mean((ptr_std_forecast_multi_2-ptr_prox4Cf).^2))
              sqrt(mean((ptr_std_forecast_multi_3-ptr_prox4Cf).^2))
              sqrt(mean((ptr_std_forecast_multi_4-ptr_prox4Cf).^2))
              sqrt(mean((ptr_std_forecast_multi_5-ptr_prox4Cf).^2))]
          
% mad_multi
% Proxy1
std_mad_multi_prox1(:,5)=[mean(abs(ptr_std_forecast_multi_1-ptr_prox1Cf))
                    mean(abs(ptr_std_forecast_multi_2-ptr_prox1Cf))
                    mean(abs(ptr_std_forecast_multi_3-ptr_prox1Cf))
                    mean(abs(ptr_std_forecast_multi_4-ptr_prox1Cf))
                    mean(abs(ptr_std_forecast_multi_5-ptr_prox1Cf))]
% Proxy2
std_mad_multi_prox2(:,5)=[mean(abs(ptr_std_forecast_multi_1-ptr_prox2Cf))
                    mean(abs(ptr_std_forecast_multi_2-ptr_prox2Cf))
                    mean(abs(ptr_std_forecast_multi_3-ptr_prox2Cf))
                    mean(abs(ptr_std_forecast_multi_4-ptr_prox2Cf))
                    mean(abs(ptr_std_forecast_multi_5-ptr_prox2Cf))]
% Proxy3
std_mad_multi_prox3(:,5)=[mean(abs(ptr_std_forecast_multi_1-ptr_prox3Cf))
                    mean(abs(ptr_std_forecast_multi_2-ptr_prox3Cf))
                    mean(abs(ptr_std_forecast_multi_3-ptr_prox3Cf))
                    mean(abs(ptr_std_forecast_multi_4-ptr_prox3Cf))
                    mean(abs(ptr_std_forecast_multi_5-ptr_prox3Cf))]
% Proxy4
std_mad_multi_prox4(:,5)=[mean(abs(ptr_std_forecast_multi_1-ptr_prox4Cf))
                    mean(abs(ptr_std_forecast_multi_2-ptr_prox4Cf))
                    mean(abs(ptr_std_forecast_multi_3-ptr_prox4Cf))
                    mean(abs(ptr_std_forecast_multi_4-ptr_prox4Cf))
                    mean(abs(ptr_std_forecast_multi_5-ptr_prox4Cf))]

%% Final Sort
% rtn
 [rtn_rmse_multi_ranked_val,rtn_rmse_multi_rank]=sort(rtn_rmse_multi,1);
 [rtn_mad_multi_ranked_val,rtn_mad_multi_rank]=sort(rtn_mad_multi,1);
 % std
 % rmse_multi
 % Proxy 1
[std_prox1_rmse_multi_ranked_val,std_prox1_rtn_rmse_multi_rank]=sort(std_rmse_multi_prox1,1);
% Proxy 2
[std_prox2_rmse_multi_ranked_val,std_prox2_rtn_rmse_multi_rank]=sort(std_rmse_multi_prox2,1);
% Proxy 3
[std_prox3_rmse_multi_ranked_val,std_prox3_rtn_rmse_multi_rank]=sort(std_rmse_multi_prox3,1);
% Proxy 4
[std_prox4_rmse_multi_ranked_val,std_prox4_rtn_rmse_multi_rank]=sort(std_rmse_multi_prox4,1);
% mad_multi
 % Proxy 1
[std_prox1_mad_multi_ranked_val,std_prox1_rtn_mad_multi_rank]=sort(std_mad_multi_prox1,1);
% Proxy 2
[std_prox2_mad_multi_ranked_val,std_prox2_rtn_mad_multi_rank]=sort(std_mad_multi_prox2,1);
% Proxy 3
[std_prox3_mad_multi_ranked_val,std_prox3_rtn_mad_multi_rank]=sort(std_mad_multi_prox3,1);
% Proxy 4
[std_prox4_mad_multi_ranked_val,std_prox4_rtn_mad_multi_rank]=sort(std_mad_multi_prox4,1);
