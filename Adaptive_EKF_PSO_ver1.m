clear 
close all
% Load ECG data from MAT file
[file, path] = uigetfile('*.mat','Select ECG mat file');
data = load(file);


fs = data.fs;                % Sampling frequency of ECG signal (Hz)

x = data.x;               % Extract stored signal matrix
ecg = x(1,:);             % Use first channel as ECG signal
length_sig = length(ecg); % Total number of ECG samples

ecg_bins = round(fs/2);           % Number of phase bins for ECG mean calculation

dt = 1/fs; % time step

%% adding white gaussian Noise
SNR = 6;   % you can change the Noise SNR level
x_noisy = zeros(3,length(x(1,:)));
x_noisy(1,:) = awgn(x(1,:),SNR,'measured');
%% -------- R-peak detection using Pan–Tompkins algorithm
[qrs_positions] = pantompkins_qrs(x_noisy(1,:),fs);
figure(1),plot(x_noisy(1,:),'b'),hold on,plot(qrs_positions,x_noisy(1,qrs_positions),'*r'),hold off
legend({'Noisy ECG Signal','R Peaks'})
title([file '   at SNR = ' num2str(SNR)])
axis('tight')
%% -------- Phase calculation
% Linear phase based on RR intervals
[Linearphase,~] = calculate_linear_phase_ver2(qrs_positions,length_sig,fs);
x_noisy(2,:) = Linearphase; 



[ECGsd,ECGmean,meanphase] = ecgsd_extractor_ver1(x_noisy(1,:),Linearphase,ecg_bins);

% further smoothing of ECG mean using wavelet
ECGmean = wdenoise(ECGmean,5,Wavelet="bior4.4",DenoisingMethod="BlockJS");



%% -------- ECG parameter extraction using Gaussian mixture model

MaxNumGaussian = 50;   % Maximum number of Strongest Gaussian components

% ========================building new myfun based on L Gaussians
L_num_of_Gaussian_kernels = 50;
ecg_mean_temp = 0;
ai = [];
bi = [];
tetai  = [];
for i=1:L_num_of_Gaussian_kernels
% disp(num2str(i))
ecg_mean_temp1 = ECGmean - ecg_mean_temp;
lb = [-1.5*max(ecg_mean_temp1).*ones(1,1)   0.000001*ones(1,1)   (-pi+.014)*ones(1,1)  ];
ub = [(1.5*max(ecg_mean_temp1)).*ones(1,1)  5*ones(1,1)  (pi-.014)*ones(1,1)  ];  
myfun1 = @(params)  norm(ecg_mean_temp1'-sum((repmat(params(1:1),ecg_bins,1).*exp(-(rem(repmat(meanphase,1,1)'-repmat(params(3),ecg_bins,1)+pi,2*pi)-pi) .^2 ./ (2*(repmat(params(2),ecg_bins,1)) .^ 2))),2));


% options = optimoptions('particleswarm','SwarmSize',30,'HybridFcn',@fmincon,'MaxIter',1000);
options = optimoptions('particleswarm','SwarmSize',50,'MaxIter',100,'Display','off');

OptimumParams = particleswarm(myfun1,3*1,lb,ub,options);

% L = (length(OptimumParams)/3);

ai_1 = OptimumParams(1);
bi_1 = OptimumParams(2);
tetai_1 = OptimumParams(3);
ai = [ai ai_1];
bi = [bi bi_1];
tetai  = [tetai tetai_1];
dtetai_1 = rem(meanphase - tetai_1 + pi,2*pi)-pi;
ecg_mean_temp = ecg_mean_temp + ai_1 .* exp(-dtetai_1 .^2 ./ (2*bi_1 .^ 2));
figure(41),plot(ecg_mean_temp,'b'),hold on,plot(ECGmean,'r')
legend({'Synthetic ECG','ECG Mean'}),hold off
title([num2str(i) 'th' '  Gaussian Kernel found'])
axis tight
% pause(3)
end

%% selection of the Strongest Peaks
[~,indx_strongest_peaks] = sort(abs(ai),'descend');

ai = ai(indx_strongest_peaks(1:MaxNumGaussian));
bi = bi(indx_strongest_peaks(1:MaxNumGaussian));
tetai = tetai(indx_strongest_peaks(1:MaxNumGaussian));


Alpha_i = ai;   % Gaussian amplitudes
Beta_i  = bi;   % Gaussian widths
Theta_i = tetai;   % Gaussian centers


% Sorting of parameters from based on tetai from -pi to pi

[Theta_i,idx] = sort(Theta_i,'ascend');
Alpha_i = Alpha_i(idx);
Beta_i = Beta_i(idx);
OptimumParams = [Alpha_i Beta_i Theta_i];
params = OptimumParams;
size_params = length(params);


%%  ANGULAR FREQUENCY MEASUREMENT

ind=1*qrs_positions;
ind2 = ind-[0 ind(1:end-1)];
RR = mean(ind2(1,2:end)); % mean of RR Intervals
w=2*pi*fs/RR; % angular frequency
RR_var=std(2*pi*fs./ind2(1,2:end));% standard deviation of RR Intervals







RR = mean(diff(ind(2:end-1)));

stepteta=2*pi/RR;
w_1 = fs*stepteta;
for j=ind(1,1):-1:1
    
     x_noisy(3,j) = w_1;

end   
for i=1:length(ind)-1



bins = ind(1,i+1)-ind(1,i);

 stepteta = 2*pi/(bins);
w_1 = fs*stepteta;
for j=ind(1,i)+1:ind(1,i+1)
 
    x_noisy(3,j) = w_1;

   
end

end
stepteta=2*pi/RR;
teta= 0;
w_1 = fs*stepteta;

for j=min(ind(end,end)+1,size(x_noisy,2)):size(x_noisy,2)

        x_noisy(3,j) = w_1;

 
end



%% EKF3 Initializations
ai = params(1,1:size_params/3);
bi = params(1,size_params/3+1:2*size_params/3);
tetai = params(1,2*size_params/3+1:size_params);
 

RR_var = var(x_noisy(3,:));

R = diag([ 0.1*mean(ECGsd)^2 (.001) RR_var ]);
Q = diag( [   .01*R(1,1) (.1) 0.01*RR_var (.05.*ai.*ones(1,size_params/3)).^2 (.05*bi.*ones(1,size_params/3)).^2 (.005*tetai.*ones(1,size_params/3)).^2  ] );



y = [x_noisy(1,:);x_noisy(2,:);x_noisy(3,:)]; % noisy measurements 
 


% creating output vectors Xeks Xekf Xpred
Xekf_update = zeros(3+0,size(x_noisy,2));
Xekf_pred = zeros(size(Xekf_update));




H = eye(3);
Ak = eye(0+3);
Fk = zeros(3,size_params+3);
Fk(1:3,1:3)=eye(3);
Ak(2,3)=dt;

%% new lines
Ak = repmat(Ak,[1,1,size(x_noisy,2)]);
Fk = repmat(Fk,[1,1,size(x_noisy,2)]);


%%
X0=[x_noisy(1,1);x_noisy(2,1);w]; 

P0 = 0.1*eye(3);
P0(1,1) = 10;
Pp0  = repmat(P0,[1,1,size(x_noisy,2)]);
Pp =Pp0;

ai = params(1,1:size_params/3);
bi = params(1,size_params/3+1:2*size_params/3);
tetai = params(1,2*size_params/3+1:size_params);
 






%% Covariance Adaptation matrix parameters
Memory_R = [];
window_size = round(fs/3);
forgetting_factor = 0.99;
%% EKF filtering
for i=1:size(Xekf_update,2)



     Pp0(:,:,i) = P0;
     Xekf_pred(:,i)= X0;
     
     %%%% updating%%%%%%%%%%
     
     

    
    
   K=P0*H'/(H*P0*H'+R);
   P=P0-K*H*P0+eps; %% updated covariance
   P=(P+P')/2;

   Xekf_update(:,i)= X0+K*(y(:,i)-H*X0);
    
   
    %%%%%%%%%%%%% Jacobian of transition matrix (A) and noise matrix(F)%%%%%
   %%% /////////Jacobian of transition matrix (A)\\\\\\\\\\\
   w = Xekf_update(3,i)'; 
   Memory_R = [Memory_R y(:,i)-H*X0];
   error = y(:,i)-H*X0;
   
   %% Adaptation of measurement and state Covariance matrices
   if i>(window_size+length(X0)+1)
       R(1,1) = (1-forgetting_factor)*mean(Memory_R(1,end-window_size+1:end-1).^2)+forgetting_factor*R(1,1);

       temp = K*(error*error')*K';
       Q(1,1)= (1-forgetting_factor)*Q(1,1)+(forgetting_factor)*temp(1,1);  % it is working
       Memory_R(:,1) = [];

   end


   dtetai = rem(Xekf_update(2,i)-tetai+pi,2*pi)-pi;
   
   
   %% new lines
   Ak(1,2,i) = -dt*sum( w*ai(1:end)./(bi(1:end).^2).*(1 - dtetai(1:end).^2./bi(1:end).^2).*exp(-dtetai(1:end).^2./(2*bi(1:end).^2)) ) ; %%df1/teta
    

   Aw = - dt*ai.*dtetai./(bi.^2).*exp(-dtetai.^2./(2*bi.^2));
   Ak(1,3,i) = sum(Aw(1:end));%%df1/w
   
   
   
   
  
   
   
   
   %%%%%% derivation with respect to ai
   Fa =-dt*w./(bi.^2).*dtetai .* exp(-dtetai.^2./(2*bi.^2));
   Fk(1,4:size_params/3+3,i) = Fa(1:end) ;%% dF1/ai
   
   
   %%%%%% derivation with respect to bi
   Fb = 2*dt.*ai.*w.*dtetai./bi.^3.*(1 - dtetai.^2./(2*bi.^2)).*exp(-dtetai.^2./(2*bi.^2));
   
   Fk(1,size_params/3+4:2*size_params/3+3,i)= Fb(1:end);%% dF1/bi

   
  %%%%%% derivation with respect to tetai
   Ft = dt*w*ai./(bi.^2).*exp(-dtetai.^2./(2*bi.^2)) .* (1 - dtetai.^2./bi.^2);
   Fk(1,2*size_params/3+4:size_params+3,i) = Ft(1:end);  %% dF1/tetai
   
  
   
   
   
   

   
   
   
   

   
   
   %%%%%%%%% prediction%%%%%%%%%%%%%%
   X0(1,1) = Xekf_update(1,i)-dt*sum(w*ai./(bi.^2).*dtetai.*exp(-dtetai.^2./(2*bi.^2))); 
   X0(2,1) = rem(Xekf_update(2,i)+w*dt,2*pi);
   X0(3,1)=w;
   
   
   P0 = Ak(:,:,i)*P*Ak(:,:,i)'+Fk(:,:,i)*Q*Fk(:,:,i)';
   %%%% storing updated covariance matrix
   Pp(:,:,i) = P;
  

end





%% smoothing
Peks =  Pp;
Xeks = Xekf_update;
for i=size(Xekf_update,2)-1:-1:1

    
    
   S = Pp(:,:,i) * Ak(:,:,i)' /(Pp0(:,:,i+1));
   Xeks(:,i) = Xekf_update(:,i) + S * (Xeks(:,i+1) - Xekf_pred(:,i+1));
   Peks (:,:,i) = Peks(:,:,i) - S * (Pp0(:,:,i+1) - Peks(:,:,i+1)) * S';



end




AEKF_SNR = 10*log10(mean((x(1,:)-x_noisy(1,:)).^2)/mean((x(1,:)-Xekf_update(1,:)).^2))
AEKS_SNR = 10*log10(mean((x(1,:)-x_noisy(1,:)).^2)/mean((x(1,:)-Xeks(1,:)).^2))


figure(2),
subplot(3,1,1)
        plot(1:length(x),x(1,:),'k')
                legend({'Original'})
                        axis('tight')
        title([file])
        subplot(3,1,3)
        plot(1:length(x),Xekf_update(1,:),'r')
                legend({'Adaptive EKF'})
                        axis('tight')

        subplot(3,1,2),plot(1:length(x),x_noisy,'b')
        legend({'Noisy'})
        axis('tight')


figure(3),
subplot(3,1,1)
        plot(1:length(x),x(1,:),'k')
                legend({'Original'})
                        axis('tight')
        title([file])
        subplot(3,1,3)
        plot(1:length(x),Xeks(1,:),'r')
                legend({'Adaptive EKS'})
                        axis('tight')

        subplot(3,1,2),plot(1:length(x),x_noisy,'b')
        legend({'Noisy'})
        axis('tight')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
function [Phase,Omega] = calculate_linear_phase_ver2(locs,length_sig,fs)

% locs       : indices of detected R-peaks
% length_sig : total number of ECG samples
% fs         : sampling frequency

ind = locs(:)';                % Convert R‑peak indices to row vector

Phase = zeros(1,length_sig);   % Phase of each ECG sample
Omega = zeros(1,length_sig);   % Instantaneous angular frequency

RR = mean(diff(ind));          % Mean RR interval (samples)

%% -------- Phase before the first R‑peak

stepTheta = 2*pi/RR;           % Average phase increment per sample
omega_val = fs*stepTheta;      % Instantaneous angular frequency

theta = 0;                     % Initialize phase

for j = ind(1)-1:-1:1          % Move backward from first R‑peak
    theta = theta - stepTheta; % Decrease phase
    theta = mod(theta+pi,2*pi)-pi; % Wrap phase into [-pi , pi]

    Phase(j) = theta;          % Store phase
    Omega(j) = omega_val;      % Store frequency
end

%% -------- Phase between consecutive R‑peaks

for k = 1:length(ind)-1

    bins = ind(k+1)-ind(k);    % Number of samples between R-peaks

    stepTheta = 2*pi/bins;     % Phase increment so phase spans one cycle
    omega_val = fs*stepTheta;  % Corresponding angular frequency

    theta = 0;
    Phase(ind(k)) = 0;         % Define phase at R‑peak as zero

    for j = ind(k)+1 : ind(k+1)-1
        theta = theta + stepTheta; % Linear phase progression
        if theta>pi
            theta = -pi;
        end
        Phase(j) = theta;
        Omega(j) = omega_val;
    end

    Phase(ind(k+1)) = 0;       % Next R‑peak also set to zero phase
end

%% -------- Phase after the last R‑peak

stepTheta = 2*pi/RR;           % Use mean RR again
omega_val = fs*stepTheta;

theta = 0;

for j = ind(end)+1:length_sig
    theta = theta + stepTheta; % Continue phase linearly
    theta = mod(theta+pi,2*pi)-pi;

    Phase(j) = theta;
    Omega(j) = omega_val;
end

end



function [ecgsd,ecg_mean,phase_mean] = ecgsd_extractor_ver1(ecg,phase,bins)

x1 = ecg;                        % ECG signal
meanPhase = zeros(1,bins);       % Mean phase per bin
ECGmean = zeros(1,bins);         % Mean ECG per bin
ECGsd = zeros(1,bins);           % ECG standard deviation per bin

% Handle wrap-around phase bin near -pi / +pi
I = find( phase >= (pi-pi/bins) | phase < (-pi+pi/bins) );

if(~isempty(I))
    meanPhase(1) = -pi;
    ECGmean(1) = mean(x1(I));
    ECGsd(1) = std(x1(I));
else
    ECGsd(1) = -1;               % Mark empty bins
end

% Loop over phase bins
for i = 1 : bins-1
    I = find( phase >= 2*pi*(i-0.5)/bins - pi & ...
              phase <  2*pi*(i+0.5)/bins - pi );

    if(~isempty(I))
        meanPhase(i+1) = mean(phase(I));
        ECGmean(i+1) = mean(x1(I));
        ECGsd(i+1) = std(x1(I));
    else
        ECGsd(i+1) = -1;
    end
end

% Interpolate missing bins
K = find(ECGsd==-1);

for i = 1:length(K)
    switch K(i)
        case 1
            meanPhase(1) = -pi;
            ECGmean(1) = ECGmean(2);
            ECGsd(1) = ECGsd(2);
        case bins
            meanPhase(bins) = pi;
            ECGmean(bins) = ECGmean(bins-1);
            ECGsd(bins) = ECGsd(bins-1);
        otherwise
            meanPhase(K(i)) = mean(meanPhase([K(i)-1 K(i)+1]));
            ECGmean(K(i))   = mean(ECGmean([K(i)-1 K(i)+1]));
            ECGsd(K(i))     = mean(ECGsd([K(i)-1 K(i)+1]));
    end
end

phase_mean = meanPhase;
ecg_mean   = ECGmean;
ecgsd      = ECGsd;

end