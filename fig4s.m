%Physics567 final presentation fig4 - soma f/I curve
% run for data
clear;clc;
global RT RS RD CS CD Vrd Vrs theta Eca reset

%% Parameters
RT = 65; % in Mohm %model distal dendrite (original 65)
RS = 50; % in Mohm
RD = 43; % in Mohm
CS = 0.26; % in nF
CD = 0.12; % in nF
Vrs = -70; % in mV
Vrd = -60; % in mV
theta = -47; %in mV (soma threshold)
reset = -52; % in mV (soma reset)

% potassium current
gAHP = 4E-6 ; % in mS
EK = -90 ; % in mV
tauK = 80 ;% in ms

%calcium current
gCa = 70E-6 ;% in mS
taum = 15 ;% in ms
tauh = 80 ;% in ms

% time vector in ms
t1 = 0 ;
t2 = 48000; 
dt = 0.1;
tvec = t1:dt:t2;

%vectors (vs,vd,m,h,Iahp,ICa)
vs = zeros(size(tvec));
vd = zeros(size(tvec));
vs(1) = Vrs;
vd(1) = Vrd;
m = ones(size(tvec));
h = zeros(size(tvec));
Iahp = zeros(size(tvec));
ICa = zeros(size(tvec));

Eca = 137;  % Equilibrium potential for Ca 
            % [https://www.physiologyweb.com/lecture_notes/resting_membrane_potential/resting_membrane_potential_nernst_equilibrium_potential.html]

%% Noisy input current

sigma = 0.3; %nA
tau = 3; %ms

% mus = [0, 0.03, 0.06, 0.09, 0.12, 0.15,...
%     0.18, 0.21, 0.24, 0.27, 0.30, 0.33]; %in nA
% test soma threshold
mus = 0:0.05:0.55;%step 100pA
% injd
% mus = 1:0.06:2.3;
            
%% iterations
cvs=[];
num_iters = 1;
freqs_sum = zeros(1, length(mus));
for i = 1:num_iters
    j=0;
    %generate noisy current imput
    currents1 = [];
    muvec = [];
    for mu = mus
        current = [];
        current(1) = mu;
        for i = 2:(length(tvec) / length(mus))
            current(i) = current(i-1) + ((mu - current(i-1)) / tau)*dt + ...
                sigma*normrnd(0,1)*sqrt((2*dt) / tau);
           
        end
        muv = ones(size(current)).*mu;
        currents1 = horzcat(currents1, current);
        muvec = horzcat(muvec,muv);
    end
    
    for i = length(currents1):length(tvec)
        currents1(i) = mu(end);
    end
    movave1 = movmean(currents1, 1000);

%     currents2 = [];
%     for mu = mus
%         current = [];
%         current(1) = mu;
%         for i = 2:(length(tvec) / length(mus))
%             current(i) = current(i-1) + ((mu - current(i-1)) / tau)*dt + ...
%                 sigma*normrnd(0,1)*sqrt((2*dt) / tau);
%         end
%         currents2 = horzcat(currents2, current);
%     end
%     for i = length(currents2):length(tvec)
%         currents2(i) = mu(end);
%     end
%     movave2 = movmean(currents2, 1000);
    
    Iinjs = currents1;
    %Iinjd = currents2;
    Iinjd = zeros(1, length(currents1));
    
    for i = 1:length(tvec)
        if vs(i) == reset
            j = j+1;
            tsum = j*dt;
            vd(i+1) = vd(i) + 10; %elevated by 10 mV
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i)); % 10^-3V * 10^-6 A/V = 10^-9 amps
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT + Iahp(i) + Iinjs(i));  
        elseif vs(i) < theta && j==0
            tsum=0;
            Iahp(i) = 0; %gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));    
       elseif vs(i) < theta && j>0
            j =j+1;
            tsum = j*dt;
            Iahp(i) =gAHP*(EK - vs(i))*exp(-tsum/tauK) ; 
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));  
        elseif vs(i) == 10
            j = j+1;
            tsum = j*dt;
            vs(i+1) = reset;
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));
        else
            j=1;
            tsum = j*dt;
            vs(i+1) = 10; %for 1ms
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vd(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vd(i)));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));   
        end
    end
    window_size = floor(length(tvec) / length(mus)); %step duration of 2s
    num_spikes_s = []; %soma

    cv = [];
    loc_s = cell(size(mus));
    isi = cell(size(mus));
    for i = 0:length(mus)-1
        num_spikes_s(i+1) = nnz(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
 
        %num_spikes(i+1) = length(findpeaks(vs((i*window_size+1):((i+1)*window_size+1))))
        loc_s{i+1} = find(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
        isi{i+1} = diff(loc_s{i+1}).*dt;% in ms
        cv(i+1) = std(isi{i+1})/mean(isi{i+1});%getCV(isi{i+1});
    end
    freqs_s = ((num_spikes_s * 1000) ./ (window_size));
    cvs = vertcat(cvs,cv);
    freqs_sum = freqs_sum + freqs_s;
end
ave_freqs_s = freqs_sum ./ num_iters; %soma 
%ave_freqs1 = freqs_sum ./ num_iters; %dendrite

%% plots

%figure4A - fre to I relationship
figure(2)
%mus0 = 0:0.03:0.33;
mus_s = 0:0.05:0.55;%soma ploting
plot(mus_s, ave_freqs_s./10, 'ok','MarkerFaceColor','k')
hold on
plot(mus_d,ave_freqs_d./10,'or','MarkerFaceColor','r')
plot(mus_sd,ave_freqs_sd./10,'ok','MarkerFaceColor','r')
xlabel('\mus (nA)')
ylabel('Mean spike rate (AP/s)')
title('Current to rate transfer functions')
legend('s','d','s+d','Location','Northwest')

% hold on
% %inject current to dendrite
% plot(pamus1, ave_freqs1, '-or','MarkerFaceColor','r')
% legend('somatic injection','dendritic injection')
% xlabel('mean current (pA)')
% ylabel('mean spike rate (AP/S)')
% title('Current to rate transfer functions')
% hold off

%for dendrite use cvs_for_plotting1

cvs_for_plotting_s = nanmean(cvs, 1); %soma_cv
cvs_for_plotting_s(isnan(cvs_for_plotting_s)) = 0;

%filter out the zeros
a_cvs = nonzeros(cvs_for_plotting_s);
a_cvs = a_cvs';

s_mus = mus_s(4:end);

%figure4B
figure(4) %%CV of isi 
plot(s_mus,a_cvs,'ok','MarkerFaceColor','k')
hold on
plot(mus_sd,cvs_for_plotting_sd,'ok','MarkerFaceColor','r')
hold off
ylim([0 3])
ylabel('C.V. of ISI')
xlabel('\mus (nA)')
legend('s','s+d')
title('C.V. of mean spike rates')

% %figure2C&D - inject current to soma/dendrite alone
% %same input to soma and dendrite (1.0-1.1nA)
% 
% figure()
% plot(tvec((10*window_size+1):(11*window_size+1)),vs1((10*window_size+1):(11*window_size+1)),'r-')
% ylabel('Somatic Voltage (mV)')
% xlabel('Time(ms)')
% title('Dendritic injection')
% % figure(1)
% % tiledlayout(2,1)
% % nexttile
% % plot(tvec,vs(1:length(tvec)),'k-')
% % ylabel('Somatic Voltage (mV)')
% % title('Figure1B.Current injection at the soma')
% % nexttile
% % plot(tvec,Iinjs,'k-', tvec(1:length(muvec)), muvec,'w-','LineWidth', 1)
% % %plot(tvec,Iinjs,'k-', tvec, movave1,'b--','LineWidth', 0.8)
% % ylabel('Noisy somatic input (pA)')
% % xlabel('Time (ms)')
% 
% %figure1C - fre to I relationship
% figure(1)
% %mus0 = 0:0.03:0.33;
% mus0=0:0.06:1.2;
% mus1 =  1:0.06:2.3;
% pamus = mus0.*(10^3);
% pamus1 = mus1.*(10^3);
% plot(pamus, ave_freqs, 'ok','MarkerFaceColor','k')
% xlabel('mean current (pA)')
% ylabel('mean spike rate (AP/S)')
% hold on
% %inject current to dendrite
% plot(pamus1, ave_freqs1, 'or','MarkerFaceColor','r')
% legend('somatic injection','dendritic injection')
% hold off
% title('Current to rate transfer functions')

%figure 1D - gain of soma/dendrite
%soma
loc_soma = 0; %in um
slope0 = [];
for i = 1:length(ave_freqs)
    g0 = polyfit(pamus(i),ave_freqs(i),1);
    slope0(i) = g0(1);
end
slope0 = max(slope0);
%dendrite
loc_dendrite = 500; %in um
slope1 = [];
for i = 1:length(ave_freqs1)
    g1 = polyfit(pamus(i),ave_freqs(i),1);
    slope1(i) = g1(1);
end
slope1 = max(slope1);
g1 = polyfit(pamus1,ave_freqs1,1);
slope1 = g1(1); %find the biggest slope?
figure(3)
plot(loc_soma,slope0,'-ok','MarkerFaceColor','k')
hold on
plot(loc_dendrite,slope1,'-or','MarkerFaceColor','r')
hold off
xlim([-20 700]);xlabel('Distance from soma (um)');
ylim([0 0.12]);ylabel('Gain(AP/s/pA');
legend('somatic injection','dendritic injection')
title('Mean gain')

% plot current
% figure(1)
% tiledlayout(1,2)
% nexttile
% plot(tvec, currents1, 'c-', tvec, movave1, 'r--', 'LineWidth', 0.8)
% nexttile
% plot(tvec, currents2, 'c-', tvec, movave2, 'r--', 'LineWidth', 0.8)

% figure(2)
% plot(tvec,vs(1:length(tvec)))
% %xlim([0 100])
% title('Somatic potential')

figure(3)
plot(tvec,vd(1:length(tvec)))
%xlim([0 100])
title('Dendritic potential')

figure(4)
plot(tvec,vs(1:length(tvec)),tvec,vd(1:length(tvec)),'r')
title('Somatic and dendritic potentials')
xlim([12000 14000])

figure(5)
plot(tvec,Iahp)
title('I_{AHP}')

figure(6)
plot(mus, ave_freqs, 'r-')

figure(7)
bar(mus, ave_freqs)

%for dendrite use cvs_for_plotting1
mus_s = mus;
cvs_for_plotting = nanmean(cvs, 1);
cvs_for_plotting(isnan(cvs_for_plotting)) = 0;
figure(8) %%CV of isi 
plot(mus,cvs_for_plotting,'or','MarkerFaceColor','r')
hold on
