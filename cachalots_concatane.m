% Samuel Pinson 13/09/22

clear all 
close all



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

folder=['P083'];
% Files= dir([folder '\**\*.wav']);
Files= dir([folder '/*.wav']);
NFiles=length(Files);
% path_='data/';
% Files= dir([path_ '*.wav']);
NFiles=length(Files);


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


Nsub=2;
for subFnum=1:Nsub
for Fnum=round(NFiles*(subFnum-1)/Nsub)+1:round(NFiles*subFnum/Nsub)
    Fnum
    fname=Files(Fnum).name(1:end-4);
    path_=folder;%Files(Fnum).folder;
    str=['load ' path_ '/P0P1_' fname ];
    eval(str)

    
    if OK==1
       tP0_= [tP0_ tP0];
       tP1_= [tP1_ tP1];
       AP0_= [AP0_ AP0]; 
       AP1_= [AP1_ AP1]; 
       indP0_= [indP0_ indP0]; 
       indP1_= [indP1_ indP1]; 
       matrep_= [matrep_ matrep.']; 
       matcoh_= [matcoh_ matcoh.']; 
       mattau_coh_= [mattau_coh_ mattau_coh.']; 
       
       tau_coh_P0_= [tau_coh_P0_ tau_coh_P0]; 
       date_= [date_ date+tP0/3600/24]; 
       noise_= [noise_ noise*ones(1,Nclick)]; 
    end
    
    
end
%%
delta_click_=(date_(2:end)-date_(1:end-1))*3600*24;
silences=find(delta_click_>30);
%%

for n=1:length(silences)
    n
    indn=silences(n)+(n-1);
   tP0_= [tP0_(1:indn) tP0_(indn) tP0_(indn+1:end)];
   tP1_= [tP1_(1:indn) tP1_(indn) tP1_(indn+1:end)];
   AP0_= [AP0_(1:indn) AP0_(indn) AP0_(indn+1:end)]; 
   AP1_= [AP1_(1:indn) AP1_(indn)  AP1_(indn+1:end) ]; 
   indP0_= [indP0_(1:indn) indP0_(indn) indP0_(indn+1:end)]; 
   indP1_= [indP1_(1:indn) indP1_(indn) indP1_(indn+1:end)]; 
   matrep_= [matrep_(:,1:indn) 0*matrep_(:,indn) matrep_(:,indn+1:end)]; 
   matcoh_= [matcoh_(:,1:indn) 0*matcoh_(:,indn) matcoh_(:,indn+1:end)]; 
   mattau_coh_= [mattau_coh_(:,1:indn) 0*mattau_coh_(:,indn) mattau_coh_(:,indn+1:end)]; 

   tau_coh_P0_= [tau_coh_P0_(1:indn) tau_coh_P0_(indn) tau_coh_P0_(indn+1:end)]; 
   date_= [date_(1:indn) date_(indn)+0.2/3600/24 date_(indn+1:end)]; 
   noise_= [noise_(1:indn) noise_(indn) noise_(indn+1:end)]; 
end

%%
% figure;hplot=pcolor(date_,t_subsig,matrep_./(AP0_(ones(1,length(t_subsig)),:)));set(hplot,'linestyle','none')
figure;
subplot(2,1,1);
hplot=pcolor(date_,t_subsig(1:4:end),20*log10(matrep_(1:4:end,:)));set(hplot,'linestyle','none')
colormap('jet')
% caxis([0 1e-3])
caxis([max(20*log10(matrep_(:)))-50 max(20*log10(matrep_(:)))-20])
datetick('x', 'HHMM');
%%
subplot(2,1,2)
plot(date_,tau_coh_P0_,'.')
datetick('x', 'HHMM');
%save XXXX t_subsig subsiglength Fs 

str=['save repmatP083_' num2str(subFnum) '_' num2str(Nsub) ' Fs date_ tau_coh_P0_ matrep_ indP0_ indP1_ tP0_ tP1_'];

eval(str)
end

%%% save repmatP083_2_2 Fs date_ tau_coh_P0_ matrep_ indP0_ indP1_ tP0_ tP1_
