%% experiment with RT - fig2C&D run for dendrite injection
%RT will affect soma response to somatic injections
clear;clc;
global RT RS RD CS CD Vrd Vrs theta Eca reset

%% Parameters
RT = 15; % in Mohm %model distal dendrite (original 65)
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
t2 = 42000;
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
 %mus = 0:0.06:1.2;%step 100pA
% injd
 mus = 0.6:0.06:1.8;
 %mus1=mus;           
%% iterations
cvs1=[];
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
%     
    Iinjd = currents1;
    %Iinjd = currents2;
    Iinjs = zeros(1, length(currents1));
    
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
    num_spikes = [];
    cv = [];
    loc = cell(size(mus));
    isi = cell(size(mus));
    for i = 0:length(mus)-1
        num_spikes(i+1) = nnz(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
        %num_spikes(i+1) = length(findpeaks(vs((i*window_size+1):((i+1)*window_size+1))))
        loc{i+1} = find(vs((i*window_size+1):((i+1)*window_size+1)) == 10);
        isi{i+1} = diff(loc{i+1}).*dt;% in ms
        cv(i+1) = std(isi{i+1})/mean(isi{i+1});%getCV(isi{i+1});
    end
    freqs = ((num_spikes * 1000) ./ (window_size));
    cvs1 = vertcat(cvs1,cv);
    freqs_sum = freqs_sum + freqs;
end
%ave_freqs = freqs_sum ./ num_iters;  
ave_freqs = freqs_sum ./ num_iters; %dendrite

vs1 = vs;

%% plot soma voltage 2D
figure(2)
%mu = 0.9

vs0=vs1((9*window_size):(9*window_size)+2*10^3);
tvec1 = tvec((9*window_size):(9*window_size)+2*10^3);
plot(tvec1,vs0,'r-')
hold on
plot(tvec4,vs4,'b-')
plot([18020; 18040], [-65; -65], '-k',  [18020; 18020], [-65; -55], '-k', 'LineWidth', 2)
text(18032,-68, '20 ms', 'HorizontalAlignment','right')
text(18012,-60, '10 mV', 'HorizontalAlignment','center')
set(gca, 'Visible', 'off')
set(gca,'fontsize',15)
legend('Dendritic injection','Location','Southeast')
hold off

%zoom in dendritic injection

figure(4)
tvec4 = tvec((9*window_size)+1*10^3:(9*window_size)+1.2*10^3);
vs4=vs1((9*window_size)+1*10^3:(9*window_size)+1.2*10^3);
plot(tvec4,vs4,'r-')
hold on 
plot([18100; 18105], [-65; -65], '-k','LineWidth', 2)
text(18103, -69,'5 ms','HorizontalAlignment','right')
set(gca, 'Visible', 'off')
set(gca,'fontsize',15)
legend('Zoom-in somatic injection','Location','Southeast')
hold off

cvs_for_plotting1 = nanmean(cvs1, 1);
cvs_for_plotting1(isnan(cvs_for_plotting1)) = 0;
figure(8) %%CV of isi 
plot(mus1,cvs_for_plotting1,'or','MarkerFaceColor','r')
hold on


