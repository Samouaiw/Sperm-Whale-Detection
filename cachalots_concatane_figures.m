% Samuel Pinson 13/09/22

clear all
close all

load repmatP092_1_2 
date=date_; 
tau_coh_P0=tau_coh_P0_ ;
matrep=matrep_ ;
indP0=indP0_ ;
indP1=indP1_ ;
tP0=tP0_ ;
tP1=tP1_;

clear date_ tau_coh_P0_ matrep_ indP0_ indP1_ tP0_ tP1_


load repmatP092_2_2 

date=[date date_]; 
tau_coh_P0=[tau_coh_P0 tau_coh_P0_] ;
matrep=[matrep matrep_] ;
indP0=[indP0 indP0_] ;
indP1=[indP1 indP1_] ;
tP0=[tP0 tP0_] ;
tP1=[tP1 tP1_];

clear date_ tau_coh_P0_ matrep_ indP0_ indP1_ tP0_ tP1_

%approximative a priori:
t_P0_P1_approx=10e-3; % s
t_P_approx=2e-3; % s
T_click=30e-3; % s
%%%%%%%%%%%%%%%%%%%%%%%%
Fs=39062;
Ndecimate=16;%
subsiglength=1; % s
t_subsig=([1:Ndecimate:round(subsiglength*Fs)]-round(T_click*Fs))/Fs;
Nt_s=length(t_subsig);
% %%
% figure(1)
% subplot(2,1,1)
% plot(date,tau_coh_P0*1000,'.')
% ylabel('\tau_{coh} (ms)')
% subplot(2,1,2)
% hplot=pcolor(date,t_subsig(1:4:round(Nt_s/2)),20*log10(matrep(1:4:round(Nt_s/2),1:end)/max(matrep(:))));set(hplot,'linestyle','none')
% colormap('jet')
% caxis([-50 -20])
% ylim([min(t_subsig) 0.3])
% % xlim([date(1) date(end)])
% % datetick('x', 'HHMM','keeplimits');

% duree=90;%minutes
% overlap=30;%minutes
% N=round((date(end)-date(1))/((duree-overlap)/24/60));
% for n=1:N
%     date1=date(1)+n*(duree-overlap)/24/60;
%     date2=date(1)+n*(duree-overlap)/24/60+duree/24/60;
%     subplot(2,1,1)
%     xlim([date1 date2])
%     datetick('x', 'HHMM','keeplimits');
%     subplot(2,1,2)
%     xlim([date1 date2])
%     datetick('x', 'HHMM','keeplimits');
%     pause
% end





%%
date1_d=datenum(2017,1,21,04,20,00);
date2_d=datenum(2017,1,21,05,00,00);
date1_b=datenum(2017,1,21,01,00,00);
date2_b=datenum(2017,1,21,01,40,00);
date1_w=datenum(2017,1,21,10,40,00);
date2_w=datenum(2017,1,21,11,00,00);
% xlim([date1 date2])
% datetick('x', 'HHMM','keeplimits');


[AA,inddate1_d]=min(abs(date1_d-date));
[AA,inddate2_d]=min(abs(date2_d-date));
[AA,inddate1_b]=min(abs(date1_b-date));
[AA,inddate2_b]=min(abs(date2_b-date));
[AA,inddate1_w]=min(abs(date1_w-date));
[AA,inddate2_w]=min(abs(date2_w-date));

% %%
% figure(2)
% subplot(3,1,1)
% plot(date,tau_coh_P0*1000,'.')
% ylabel('\tau_{coh} (ms)')
% xlim([min(date),max(date)])
% datetick('x', 'HHMM','keeplimits');
% title('(a)')
% 
% subplot(3,1,2)
% hplot=pcolor(date,t_subsig(1:4:round(end/2)),20*log10(matrep(1:4:round(end/2),:)/max(matrep(:))));set(hplot,'linestyle','none')
% colormap('jet')
% colorbar
% ylabel('Time (s)')
% caxis([-50 -20])
% xlim([min(date) max(date)])
% ylim([min(t_subsig) 0.25])
% datetick('x', 'HHMM','keeplimits');
% title('(b)')
% 
% subplot(3,3,7)
% hplot=pcolor(date(inddate1_d:inddate2_d),t_subsig,20*log10(matrep(:,inddate1_d:inddate2_d)/max(matrep(:))));set(hplot,'linestyle','none')
% colormap('jet')
% ylabel('Time (s)')
% caxis([-50 -20])
% ylim([min(t_subsig) 0.1])
% xlim([date1_d date2_d])
% datetick('x', 'HHMM','keeplimits');
% title('(c)')
% 
% subplot(3,3,8)
% hplot=pcolor(date(inddate1_b:inddate2_b),t_subsig,20*log10(matrep(:,inddate1_b:inddate2_b)/max(matrep(:))));set(hplot,'linestyle','none')
% colormap('jet')
% %xlabel('Time (s)')
% caxis([-50 -20])
% ylim([min(t_subsig) 0.25])
% xlim([date1_b date2_b])
% datetick('x', 'HHMM','keeplimits');
% title('(d)')
% 
% subplot(3,3,9)
% hplot=pcolor(date(inddate1_w:inddate2_w),t_subsig,20*log10(matrep(:,inddate1_w:inddate2_w)/max(matrep(:))));set(hplot,'linestyle','none')
% colormap('jet')
% %xlabel('Time (s)')
% caxis([-50 -20])
% ylim([min(t_subsig) 0.1])
% xlim([date1_w date2_w])
% datetick('x', 'HHMM','keeplimits');
% title('(e)')

%%
figure(3)
subplot(3,1,1)
plot(date,tau_coh_P0*1000,'.')
ylabel('\tau_{coh} (ms)')
xlim([min(date),max(date)])
datetick('x', 'HHMM','keeplimits');


subplot(3,1,2)
hplot=pcolor(date,t_subsig(1:4:round(end/2)),20*log10(matrep(1:4:round(end/2),:)/max(matrep(:))));set(hplot,'linestyle','none')
colormap('jet')
colorbar
ylabel('Time (s)')
caxis([-50 -20])
xlim([min(date) max(date)])
ylim([min(t_subsig) 0.25])
datetick('x', 'HHMM','keeplimits');


subplot(3,3,7)
hplot=pcolor(date(inddate1_d:inddate2_d),t_subsig,20*log10(matrep(:,inddate1_d:inddate2_d)/max(matrep(:))));set(hplot,'linestyle','none')
colormap('jet')
ylabel('Time (s)')
caxis([-50 -20])
ylim([min(t_subsig) 0.1])
xlim([date1_d date2_d])
datetick('x', 'HHMM','keeplimits');


subplot(3,3,8)
hplot=pcolor(date(inddate1_b:inddate2_b),t_subsig,20*log10(matrep(:,inddate1_b:inddate2_b)/max(matrep(:))));set(hplot,'linestyle','none')
colormap('jet')
%xlabel('Time (s)')
caxis([-50 -20])
ylim([min(t_subsig) 0.25])
xlim([date1_b date2_b])
datetick('x', 'HHMM','keeplimits');


subplot(3,3,9)
hplot=pcolor(date(inddate1_w:inddate2_w),t_subsig,20*log10(matrep(:,inddate1_w:inddate2_w)/max(matrep(:))));set(hplot,'linestyle','none')
colormap('jet')
%xlabel('Time (s)')
caxis([-50 -20])
ylim([min(t_subsig) 0.1])
xlim([date1_w date2_w])
datetick('x', 'HHMM','keeplimits');

