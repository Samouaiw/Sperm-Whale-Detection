% Samuel Pinson 13/09/22


clear all 
close all

max_size_fig=1; % max size figure in GigaBytes
max_size_fig=max_size_fig/8; % max size figure in double precision value

max_t_silence=30;

cw=1460;
dh1h2=4.6;
bathy=715;
prof=200;
fl=3000;
fh=12000;


%approximative a priori:
t_P0_P1_approx=10e-3; % s
t_P_approx=2e-3; % s
deltat_h1h2_max=dh1h2/cw; % s
deltat_echo_max=(prof+dh1h2)*2/cw; % s
deltat_click=0.5; % s
T_click=30e-3; % s
%%%%%%%%%%%%%%%%%%%%%%%%

folder=['P092'];
% Files= dir([folder '/**/*.wav']);
Files= dir([folder '/*.wav']);
NFiles=length(Files);
% path_='data/';
% Files= dir([path_ '*.wav']);
NFiles=length(Files);

[date_r,depth_r]=read_recorder_data(folder);


tP0_= [];
tP1_= [];
AP0_= []; 
AP1_= []; 
indP0_= []; 
indP1_= []; 
matrep_= []; 
matcoh_= []; 
mattau_coh_= []; 
tau_coh_P0_= []; 
date_= []; 
noise_= []; 

subFnum=3;
Nsub=4;

for Fnum=round(NFiles*(subFnum-1)/Nsub)+1:round(NFiles*subFnum/Nsub)
    Fnum
    fname=Files(Fnum).name(1:end-4);
    path_=Files(Fnum).folder;
    str=['load ' path_ '/P0P1_' fname ];
    eval(str)

    
    if OK==1

       Nclick=length(matrep(:,1));
       tP0_= [tP0_ tP0];
       tP1_= [tP1_ tP1];
       AP0_= [AP0_ AP0]; 
       AP1_= [AP1_ AP1]; 
       indP0_= [indP0_ indP0]; 
       indP1_= [indP1_ indP1]; 
       matrep_= [matrep_ matrep.']; 
       matcoh_= [matcoh_ matcoh.']; 
%        mattau_coh_= [mattau_coh_ mattau_coh.']; 
       
       tau_coh_P0_= [tau_coh_P0_ tau_coh_P0]; 
       date_= [date_ date+tP0/3600/24]; 
       noise_= [noise_ noise*ones(1,Nclick)]; 

    end
    
    
end
%%
delta_click_=(date_(2:end)-date_(1:end-1))*3600*24;
silences=find(delta_click_>max_t_silence);
Nsil=length(silences);
Nsigs=length(tP0_);
%%
   tP0_=[tP0_ zeros(1,Nsil)];
   tP1_=[tP1_ zeros(1,Nsil)];
   AP0_=[AP0_ zeros(1,Nsil)];
   AP1_=[AP1_ zeros(1,Nsil)];
   indP0_=[indP0_ zeros(1,Nsil)];
   indP1_=[indP1_ zeros(1,Nsil)];
   matrep_=[matrep_ zeros(length(t_subsig),Nsil)];
   matcoh_=[matcoh_ zeros(length(t_subsig),Nsil)];
   mattau_coh_=[mattau_coh_ zeros(length(t_subsig),Nsil)];

   tau_coh_P0_=[tau_coh_P0_ zeros(1,Nsil)];
   date_=[date_ zeros(1,Nsil)];
   noise_=[noise_ zeros(1,Nsil)];
   %%
for n=1:Nsil
    n
    indn=silences(n)+(n-1);
   tP0_(1:Nsigs+n)= [tP0_(1:indn) tP0_(indn) tP0_(indn+1:end-(Nsil-(n-1)))];
   tP1_(1:Nsigs+n)= [tP1_(1:indn) tP1_(indn) tP1_(indn+1:end-(Nsil-(n-1)))];
   AP0_(1:Nsigs+n)= [AP0_(1:indn) AP0_(indn) AP0_(indn+1:end-(Nsil-(n-1)))]; 
   AP1_(1:Nsigs+n)= [AP1_(1:indn) AP1_(indn)  AP1_(indn+1:end-(Nsil-(n-1))) ]; 
   indP0_(1:Nsigs+n)= [indP0_(1:indn) indP0_(indn) indP0_(indn+1:end-(Nsil-(n-1)))]; 
   indP1_(1:Nsigs+n)= [indP1_(1:indn) indP1_(indn) indP1_(indn+1:end-(Nsil-(n-1)))]; 
   matrep_(:,1:Nsigs+n)= [matrep_(:,1:indn) 0*matrep_(:,indn) matrep_(:,indn+1:end-(Nsil-(n-1)))]; 
   matcoh_(:,1:Nsigs+n)= [matcoh_(:,1:indn) 0*matcoh_(:,indn) matcoh_(:,indn+1:end-(Nsil-(n-1)))]; 
%    mattau_coh_(:,1:Nsigs+n)= [mattau_coh_(:,1:indn) 0*mattau_coh_(:,indn) mattau_coh_(:,indn+1:end-(Nsil-(n-1)))]; 

   tau_coh_P0_(1:Nsigs+n)= [tau_coh_P0_(1:indn) tau_coh_P0_(indn) tau_coh_P0_(indn+1:end-(Nsil-(n-1)))]; 
   date_(1:Nsigs+n)= [date_(1:indn) date_(indn)+max_t_silence*(2/3)/3600/24 date_(indn+1:end-(Nsil-(n-1)))]; 
   noise_(1:Nsigs+n)= [noise_(1:indn) noise_(indn) noise_(indn+1:end-(Nsil-(n-1)))]; 
end

 

%%
% figure;hplot=pcolor(date_,t_subsig,matrep_./(AP0_(ones(1,length(t_subsig)),:)));set(hplot,'linestyle','none')
figure;
% subplot(2,1,1);
% %plot(date_r,depth_r);datetick('x', 'HHMM');
% scatter(date_,tau_coh_P0_,30,log10(max([AP1_;AP0_],[],1)),'o','filled')
% datetick('x', 'HHMM');colormap('jet');colorbar;%caxis([0 0.02])
% subplot(2,1,2)
% scatter(date_,tau_coh_P0_,30,tP1_-tP0_,'o','filled')
% datetick('x', 'HHMM');colormap('jet');colorbar;caxis([0 0.012])
scatter(date_,tau_coh_P0_,5,'ok','filled')
datetick('x', 'HHMM');


% subplot(3,1,3)
% plot(date_,tP1_-tP0_,'.')
% datetick('x', 'HHMM');
%%
figure;
hplot=pcolor(date_,t_subsig(1:4:end),20*log10(matrep_(1:4:end,:)));set(hplot,'linestyle','none')
colormap('jet')
% caxis([0 1e-3])
caxis([max(20*log10(matrep_(:)))-50 max(20*log10(matrep_(:)))-20])
datetick('x', 'HHMM');
ylim([-0.03 0.5])
%%

%save XXXX t_subsig subsiglength Fs 