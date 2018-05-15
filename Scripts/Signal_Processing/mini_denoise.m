function [xdenoised,noisex]=mini_denoise(xd)
    dbw=1;          % Degree (Order )Wavelet
    gradacf=-10;
    acf_pre=1;
    ondeleta='sym8'; % Sort of Wavelet
    while ~or(gradacf>=0,dbw>9)
        disp(['DENOISING by wavelet level . . . . . ',num2str(dbw)])
        %         [xdenoised,~,~,~,~] = cmddenoise(xd,ondeleta,dbw);
        %          xdenoised=smooth(xd,SWS,'rloess'); ALTERNATIVE
        xdenoised=waveletdenoise(xd,ondeleta,dbw);
        noise=xd-xdenoised;
        acf=autocorr(noise,1);
        acf_act=abs(acf(end));
        gradacf=acf_act-acf_pre;
        acf_pre=acf_act;    % update ACF value
        dbw=dbw+1;          % Increase wavelet level
    end
    disp(' denoised <<<<<<<<<<<<<<<<<<<<<<<')
    dbw_opt=dbw-2; % Beacuse it detected +1 and increased it anyway +1
    %     [xdenoised,~,~,~,~] = cmddenoise(xd,ondeleta,dbw_opt);
    xdenoised=waveletdenoise(xd,ondeleta,dbw_opt);
    noisex=xd-xdenoised;
end    