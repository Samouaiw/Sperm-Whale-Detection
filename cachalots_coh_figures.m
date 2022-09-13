% Samuel Pinson 13/09/22

clear all
close all

tic

cw=1460;
dh1h2=4.6;
prof=200;
fl=3000;
fh=12000;

%approximative a priori:
t_P0_P1_approx=10e-3;
t_P_approx=2e-3;
deltat_h1h2_max=dh1h2/cw;
deltat_echo_max=(prof+dh1h2)*2/cw;
deltat_click=0.5;
T_click=30e-3;
%%%%%%%%%%%%%%%%%%%%%%%%

folder=['P081'];
% Files= dir([folder '/**/*.wav']);
Files= dir([folder '/*.wav']);
NFiles=length(Files)

for Fnum=NFiles-8
    Fnum
    path_=folder;%Files(Fnum).folder;
    %str=([path_ Files(Fnum).name]);

    %fname=['channelAB_2017-01-18_10-30-51_P81.wav'];
    %fname=['channelAB_2017-01-18_04-55-05'];
    fname=Files(Fnum).name(1:end-4);
    date_=sscanf(fname(11:end),'%4f-%2f-%2f_%2f-%2f-%2f')
    date=datenum(date_(1),date_(2),date_(3),date_(4),date_(5),date_(6))
    sig=audioread([path_ '/' fname '.wav']).';
    info=audioinfo([path_ '/' fname '.wav']);
    %sig=sig(:,1.2*1e6:round(length(sig(1,:))/10));
    Fs=info.SampleRate;
    Nt=length(sig(1,:));
    t=[0:Nt-1]/Fs;
    tshift=[-ceil(Nt/2):Nt-ceil(Nt/2)-1]/Fs;
    f=[0:Nt-1]*Fs/Nt;
    nh=2;

    %%

    wn= (2/Fs)*fl;
    [b,a]=butter(4,wn,'high');
    sigf1=filtfilt(b,a,sig(1,:));
    sigf2=filtfilt(b,a,sig(2,:));

    %%
    fc=(fh+fl)/2;
    BW=(fh-fl)/8;
    gauss_win=exp(-0.5*(f-fc).^2/BW^2);
    sigh1=(ifft(fft(sigf1).*gauss_win));
    sigh2=(ifft(fft(sigf2).*gauss_win));

    %%
    win=ifftshift(exp(-0.5*(tshift).^2/t_P_approx^2));
    win=ifft(abs(fft(win)))/sum(win); 

    %%
    Ntau_max=round(dh1h2*1.2/cw*Fs);
    tau=[-Ntau_max:Ntau_max]/Fs;
    Ntau=length(tau);

    %%
    [Nv,edges,bin]=histcounts(real(sigh1),1000);
    val=(edges(1:end-1)+edges(2:end))/2;
    sigma_=[1e-5:1e-5:1e-3].';
    Nv_=max(Nv)*exp(-0.5*val(ones(1,length(sigma_)),:).^2./sigma_(:,ones(1,length(val))).^2);
    [AA,indsigma]=min(sum(abs(Nv(ones(1,length(sigma_)),:)-Nv_).^2,2));
    gauss_n_=Nv_(indsigma,:);
    sigma_n=sigma_(indsigma);
    noise=sigma_n^2;
    %%
    figure(2)
    plot(val,(Nv))
    hold on
    plot(val,(gauss_n_),'k')
    hold off
    
    %%
    Ndecimate=16;
    sigh1_=downsample(sigh1,Ndecimate);
    sigh2_=downsample(sigh2,Ndecimate);
    win=downsample(win,Ndecimate);
    Win=fft(win);

    coh_=zeros(Ntau,length(sigh1_));
% tic
    for ntau=1:Ntau
        Fnum
        disp('time elapsed')
%         toc/60
        disp('percentage done:')
        ntau/Ntau*100
        
        sig2shift=downsample(circshift(sigh2.',round(tau(ntau)*Fs)).',Ndecimate);
        coh_(ntau,:)=ifft(fft(...
            sigh1_.*conj(sig2shift)).*Win )./(...
            sqrt(abs(ifft(fft(...
            sigh1_.*conj(sigh1_)).*Win )).*...
            abs(ifft(fft(...
            sig2shift.*conj(sig2shift)).*Win )))+noise/2);
    end

    [coh,indmaxcoh]=max(abs(coh_),[],1);
    tau_coh=tau(indmaxcoh);
    %clear coh_
    %%
%     figure(1);plot(t,abs(sigh1)/max(abs(sigh1))*3,'r')
%     hold on;
%     plot(t(1:Ndecimate:end),coh)
%     hold off
figure(1);plot(real(sig(1,:))/max(abs(sigh1)),'r')
    %%
%     ind1=3190500;
%     ind2=3193000;
    ind1=7997000;
    ind2=7999000;
    sigex=sig(1,ind1:ind2);
    sigfex=sigf1(1,ind1:ind2);
    sighex=sigh1(1,ind1:ind2);
    texshift=7.5e-3;
    tex=[0:length(sigex)-1]/Fs-texshift;
    Ntex=length(tex);
    
        
    figure(3);
    subplot(3,1,1);
    plot(tex*1000,real(sigex)/max(abs(sigex)),'k')
    xlim([min(tex*1000),max(tex*1000)])
    ylabel('Amplitude')
    title('(a)')
%     hold on;
%     plot(real(sighex)/max(abs(sighex)),'k')
%     hold off

    Nw=128;
    [S,fS,tS]=spectrogram(sigex(1,1:round(Ntex)), (1-cos(2*pi*(0:Nw)'/Nw))/2 , round(Nw*7/8) , [0:floor((Nw-1)/2)]*Fs/Nw ,Fs);
    FsS=1/mean(tS(2:end)-tS(1:end-1));

    figure(3);
    subplot(3,1,2);
    hplot=pcolor((tS-texshift)*1000,(fS)/1000,10*log10(abs(S.^2)));set(hplot,'linestyle','none')
    xlim([min(tex*1000),max(tex*1000)])
    ylabel('Frequency (kHz)')
    title('(b)')
    colormap('jet')
    colorbar
    caxis([max(10*log10(abs(S(:).^2)))-60    max(10*log10(abs(S(:).^2)))-0])

    subplot(3,1,3)
    plot(tex*1000,real(sighex)/max(abs(sighex)),'k')
    hold on
    hplot=plot(tex*1000,abs(sighex)/max(abs(sighex)),'r');set(hplot,'linewidth',2);
    xlim([min(tex*1000),max(tex*1000)])
    title('(c)')
    xlabel('Time (ms)')
    ylabel('Amplitude')
    hold off
%%    
    figure(5);
    subplot(2,1,1);
    plot(tex*1000,real(sigex)/max(abs(sigex)),'k')
    xlim([min(tex*1000),max(tex*1000)])
    ylabel('Amplitude')

    subplot(2,1,2);
    hplot=pcolor((tS-texshift)*1000,(fS)/1000,10*log10(abs(S.^2)));set(hplot,'linestyle','none')
    xlim([min(tex*1000),max(tex*1000)])
    ylabel('Frequency (kHz)')
    colormap('jet')
    colorbar
    caxis([max(10*log10(abs(S(:).^2)))-60    max(10*log10(abs(S(:).^2)))-0])
    
    figure(6);
    subplot(2,1,1);
    plot(tex*1000,real(sigex)/max(abs(sigex)),'k')
    xlim([min(tex*1000),max(tex*1000)])
    ylabel('Amplitude')

    subplot(2,1,2);
    plot(tex*1000,real(sighex)/max(abs(sighex)),'k')
    hold on
    hplot=plot(tex*1000,abs(sighex)/max(abs(sighex)),'r');set(hplot,'linewidth',2);
    xlim([min(tex*1000),max(tex*1000)])
    xlabel('Time (ms)')
    ylabel('Amplitude')
    hold off
    %%
    Nindex=10000;
    coh_ex=coh_(:,round(ind1/Ndecimate):round(ind2/Ndecimate)+Nindex);
    tau_cohex=tau_coh(1,round(ind1/Ndecimate):round(ind2/Ndecimate)+Nindex);
    cohex=coh(1,round(ind1/Ndecimate):round(ind2/Ndecimate)+Nindex);
    sigex=real(sigh1(1,ind1:ind2+Nindex*Ndecimate));
    tex_=[0:length(cohex)-1]/Fs*Ndecimate-texshift;
    tex__=[0:length(sigex)-1]/Fs-texshift;
%%
    figure(4)
    subplot(2,2,1)
    plot(tex__,sigex/max(sigex),'k')
    %xlabel('Time (s)')
    ylabel('Amplitude')
    xlim([min(tex_),max(tex_)])
    title('(a)')
    subplot(2,2,2)
    hplot=pcolor(tex_,tau*1000,abs(coh_ex));set(hplot,'linestyle','none')
    colormap('jet')
    caxis([0 1])
    xlim([min(tex_),max(tex_)])
    %xlabel('Time (s)')
    ylabel('\tau (ms)')
    title('(b)')
    colorbar
    
    subplot(2,2,3)
    plot(tex_,cohex,'k')
    %xlabel('Time (s)')
    ylabel('Coherence')
    xlim([min(tex_),max(tex_)])
    title('(c)')
    
    subplot(2,2,4)
    plot(tex_,tau_cohex*1000,'k')
    xlabel('Time (s)')
    ylabel('\tau_{coh} (ms)')
    xlim([min(tex_),max(tex_)])
    title('(d)')
%%    
    figure(7)
    subplot(3,1,1)
    plot(tex__,sigex/max(sigex),'k')
    %xlabel('Time (s)')
    ylabel('Coherence')
    xlim([min(tex_),max(tex_)])

    subplot(3,1,2)
    hplot=pcolor(tex_,tau*1000,abs(coh_ex));set(hplot,'linestyle','none')
    colormap('jet')
    caxis([0 1])
    xlim([min(tex_),max(tex_)])
    %xlabel('Time (s)')
    ylabel('\tau (ms)')

    colorbar
    

    figure(8)
    subplot(3,1,1)
    plot(tex__,sigex/max(sigex),'k')
    %xlabel('Time (s)')
    ylabel('Coherence')
    xlim([min(tex_),max(tex_)])
    
    subplot(3,1,2)
    plot(tex_,cohex,'k')
    %xlabel('Time (s)')
    ylabel('Coherence')
    xlim([min(tex_),max(tex_)])

    
    subplot(3,1,3)
    plot(tex_,tau_cohex*1000,'k')
    xlabel('Time (s)')
    ylabel('\tau_{coh} (ms)')
    xlim([min(tex_),max(tex_)])
    

    %%
%     pause(0.1)

%     str=['save ' path_ '\detect_' fname '.mat  coh tau_coh Ndecimate Nt noise date']
% %     str=['save detect_' fname ' coh tau_coh Ndecimate Nt noise date'];
%     eval(str)
end


duree1=toc;
save duree1 duree1






