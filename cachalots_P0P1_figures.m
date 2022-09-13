% Samuel Pinson 13/09/22

clear all 
close all

tic

cw=1460;
dh1h2=4.6;
prof=200;
fl=3000;
fh=12000;

seuilcoh=0.60;

%approximative a priori:
t_P0_P1_approx=10e-3; % s
t_P_approx=2e-3; % s
deltat_h1h2_max=dh1h2/cw; % s
deltat_echo_max=(prof+dh1h2)*2/cw; % s
deltat_click=0.5; % s
T_click=30e-3; % s
%%%%%%%%%%%%%%%%%%%%%%%%

coh_width_min=t_P0_P1_approx; % s
subsiglength=1; % s

%%%%%%%%%%%%%%%%%%%%%%%%

folder=['P081'];
% Files= dir([folder '/**/*.wav']);
Files= dir([folder '/*.wav']);
NFiles=length(Files)

for Fnum=NFiles-8
    Fnum
    path_=folder;%Files(Fnum).folder;
    clear tP0 tP1 AP0 AP1 indP0 indP1 tP0r tP1r AP0r AP1r indP0r indP1r tau_coh_P0
    str=([path_ '/' Files(Fnum).name]);

    fname=Files(Fnum).name(1:end-4);

    sig=audioread([path_ '/' fname '.wav']).';
    info=audioinfo([path_ '/' fname '.wav']);

    str=['load ' path_ '/detect_' fname];
    eval(str);
    
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
    win=ifftshift(exp(-0.5*(tshift).^2/(t_P_approx/2)^2));
    win=ifft(abs(fft(win)))/sum(win); 

    %%
    Ntau_max=round(dh1h2*1.2/cw*Fs);
    tau=[-Ntau_max:Ntau_max]/Fs;
    Ntau=length(tau);
    
    %%
    figure(1);;plot(t,abs(sigh1)/max(abs(sigh1))*3,'r')%max(sig(1,:))
    hold on;
    plot(t(1:Ndecimate:end),coh)
    title(fname,'interpreter','none')
    hold off
    
    %%
    sigh1w=ifft(fft(abs(sigh1)).*fft(win));
    %sigh2w=ifft(fft(sigh2).*fft(win));
    
    %%
    coh_=interp1(t(1:Ndecimate:end),coh,t);
    tau_coh_=interp1(t(1:Ndecimate:end),tau_coh,t);
    detect=zeros(size(coh_));
    detect(coh_>seuilcoh)=1;
    detect_front=zeros(size(coh_));
    detect_front(2:end)=sign(detect(2:end)-detect(1:end-1));
    ind_upfront=find(detect_front==1);
    ind_downfront=find(detect_front==-1);
    if length(ind_downfront)~=length(ind_upfront)
        if length(ind_downfront)+1==length(ind_upfront)
            ind_upfront=ind_upfront(1:end-1);
        elseif length(ind_downfront)==length(ind_upfront)+1
            ind_downfront=ind_downfront(2:end);
        else
            disp('bizarre, y a une couille dans le potage')
        end
    end
    
    coh_width=(ind_downfront-ind_upfront)/Fs;
    Ndetect=length(ind_downfront);
    m=1;
    for n=1:Ndetect
        if coh_width(n)>coh_width_min
%             [AP,indP]=max(abs(sigh1w(ind_upfront(n):ind_downfront(n))));            
            [AP,indP]=findpeaks(abs(sigh1w(ind_upfront(n):ind_downfront(n))),'SortStr','descend','NPeaks',2);
            if length(indP)>1
            indP=ind_upfront(n)-1+indP;
            %indbox=[indP-2*round(t_P0_P1_approx*Fs) indP+2*round(t_P0_P1_approx*Fs)];
            %if indbox(1)+round(subsiglength*Fs)<Nt && indbox(1)>0
                
            tau_coh_P0(m)=tau_coh_(indP(1));
            
%             sigh1_=sigh1w(indbox(1):indbox(2));
%             sigh1__=sigh1_;
%             sigh1__(2*round(t_P0_P1_approx*Fs)-round(t_P_approx*Fs):2*round(t_P0_P1_approx*Fs)+round(t_P_approx*Fs))=0;
%             [AP_,indP_]=max(abs(sigh1__));
%             indP_=indbox(1)-1+indP_;
%             if indP_>indP
%                 AP0(m)=AP;
%                 AP1(m)=AP_;
%                 tP0(m)=t(indP);
%                 tP1(m)=t(indP_);
%                 indP0(m)=indP;
%                 indP1(m)=indP_;
%             else
%                 AP0(m)=AP_;
%                 AP1(m)=AP;
%                 tP0(m)=t(indP_);
%                 tP1(m)=t(indP);
%                 indP0(m)=indP_;
%                 indP1(m)=indP;
%             end
            if indP(2)>indP(1)
                AP0(m)=AP(1);
                AP1(m)=AP(2);
                tP0(m)=t(indP(1));
                tP1(m)=t(indP(2));
                indP0(m)=indP(1);
                indP1(m)=indP(2);
            else
                AP0(m)=AP(2);
                AP1(m)=AP(1);
                tP0(m)=t(indP(2));
                tP1(m)=t(indP(1));
                indP0(m)=indP(2);
                indP1(m)=indP(1);
            end            
                      
%             indbox=[indP0(m) indP1(m)+round(deltat_echo_max*Fs)+round(T_click*Fs)];
% 
%             sigh1_=sigh1w(indbox(1):indbox(2));
%             sigh1_(1:round(2*T_click*Fs))=0;
%             [AP,indP]=max(abs(sigh1_));
%             sigh1_([[1:indP-2*round(t_P0_P1_approx*Fs)] ...
%                 [indP-round(t_P_approx*Fs):indP+round(t_P_approx*Fs)] ...
%                 [indP+2*round(t_P0_P1_approx*Fs):end]])=0;
%             [AP_,indP_]=max(abs(sigh1_));
%             
%             indP=indbox(1)-1+indP;
%             indP_=indbox(1)-1+indP_;
%             
%             if indP_>indP
%                 AP0r(m)=AP;
%                 AP1r(m)=AP_;
%                 tP0r(m)=t(indP);
%                 tP1r(m)=t(indP_);
%                 indP0r(m)=indP;
%                 indP1r(m)=indP_;
%             else
%                 AP0r(m)=AP_;
%                 AP1r(m)=AP;
%                 tP0r(m)=t(indP_);
%                 tP1r(m)=t(indP);
%                 indP0r(m)=indP_;
%                 indP1r(m)=indP;
%             end
            
            m=m+1
            %end
            end
        end
    end
    Nclick=m-1;
    if Nclick>1
%         %%        
%         figure(1);
%         plot(t,abs(sigh1))
%         hold on
%         scatter(tP0,AP0,20,'ob','filled')
%         scatter(tP1,AP1,20,'og','filled')
%         hold off
        t_subsig=([1:Ndecimate:round(subsiglength*Fs)]-round(T_click*Fs))/Fs;
        matrep=zeros(Nclick,length(t_subsig));
        matcoh=zeros(Nclick,length(t_subsig));
        mattau_coh=zeros(Nclick,length(t_subsig));
        for n=1:Nclick
            indbox(1)=indP0(n)-round(T_click*Fs);
            indbox(2)=indbox(1)-1+round(subsiglength*Fs);
            if indbox(1)>0
            if indbox(2)>=Nt
                sigh1w(end:indbox(2))=0;
                coh_(end:indbox(2))=0;
                tau_coh_(end:indbox(2))=0;
            end

            matrep(n,1:length(t_subsig))=sigh1w(indbox(1):Ndecimate:indbox(2));
            matcoh(n,1:length(t_subsig))=coh_(indbox(1):Ndecimate:indbox(2));
            mattau_coh(n,1:length(t_subsig))=tau_coh_(indbox(1):Ndecimate:indbox(2));
            
            end
        end
        Nclick=length(matrep(:,1));
        %%
        figure(2);
        subplot(2,1,1)
        h_plot=pcolor(tP0,t_subsig,20*log10(abs(matrep)).');set(h_plot,'linestyle','none')
        colormap('jet')
        caxis([ max(20*log10(abs(matrep(:))))-70 max(20*log10(abs(matrep(:))))-10])
%         hold on
%         plot(tP0,tP0r-tP0,'k')
        hold off
        subplot(2,1,2)
        plot(tP0,asin(tau_coh_P0*cw/dh1h2)*180/pi)
        
        
        %%
        
        ind1=7997000;
        ind2=7999000;
        sigex=sig(1,ind1:ind2);
        sigfex=sigf1(1,ind1:ind2);
        sighex=sigh1(1,ind1:ind2);
        sighwex=sigh1w(1,ind1:ind2);
        coh_ex=coh_(1,ind1:ind2);
        texshift=7.5e-3;
        tex=[0:length(sigex)-1]/Fs-texshift;
        Ntex=length(tex);
        
        figure(3);
        plot(tex*1000,abs(sighex)/max(abs(sighex)),'r')
        hold on
        plot(tex*1000,coh_ex,'--')
        hplot=plot(tex*1000,abs(sighwex)/max(abs(sighwex)),'k');set(hplot,'linewidth',2)
        hold off
        xlabel('Time (ms)')
        ylabel('Amplitude')
        xlim([min(tex*1000),max(tex*1000)])
        
        
        
        %%
        
%         
%         pause(1)
        %%
%         OK=1;
%         str=['save ' path_ '\P0P1_' fname '.mat tP0 tP1 AP0 AP1 indP0 indP1 matrep matcoh mattau_coh t_subsig subsiglength Fs tau_coh_P0 Nclick date noise OK'];
% 
% %         str=['save P0P1_' fname ' tP0 tP1 AP0 AP1 indP0 indP1 tP0r tP1r AP0r AP1r indP0r indP1r matrep matcoh mattau_coh t_subsig subsiglength Fs tau_coh_P0 Nclick date noise OK'];
% %         str=['save P0P1_' fname ' tP0 tP1 AP0 AP1 indP0 indP1 matrep matcoh mattau_coh t_subsig subsiglength Fs tau_coh_P0 Nclick date noise OK'];
%         eval(str)
    else
%         OK=0;
%         
%         str=['save ' path_ '\P0P1_' fname '.mat OK'];
%         eval(str)
    end
    
end

% 
% duree2=toc;
% save duree2 duree2