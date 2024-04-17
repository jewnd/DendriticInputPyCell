%Physics567 final presentation fig5 - soma+0dendrite f/I curve

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
gAHP = 4E-3 ; % in uS
EK = -90 ; % in mV
tauK = 80 ;% in ms

%calcium current (in nA)
gCa = 70E-3 ;% in uS
taum = 15 ;% in ms
tauh = 80 ;% in ms

% time vector in ms
t1 = 0 ;
t2 = 40000; 
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

sigma_s = 0.3; %nA
sigma_d = 0;
tau = 3; %ms

% mus = [0, 0.03, 0.06, 0.09, 0.12, 0.15,...
%     0.18, 0.21, 0.24, 0.27, 0.30, 0.33]; %in nA
% test soma threshold
d_inj = 0; %in nA
mus_s = 0:0.05:0.95;%step current somatic injection
mus_d = ones(1,length(mus_s)).*d_inj; %steady dendritic injection
% injd
% mus = 1:0.06:2.3;
            
%% iterations
cvs=[];
num_iters = 10;
freqs_sum = zeros(1, length(mus_s));
for i = 1:num_iters
    j=0;
    %generate noisy current imput
    currents1 = [];
    muvec = [];
    for mu = mus_s
        current = [];
        current(1) = mu;
        for i = 2:(length(tvec) / length(mus_s))
            current(i) = current(i-1) + ((mu - current(i-1)) / tau)*dt + ...
                sigma_s*normrnd(0,1)*sqrt((2*dt) / tau);
           
        end
        muv = ones(size(current)).*mu;
        currents1 = horzcat(currents1, current);
        muvec = horzcat(muvec,muv);
    end
    
    for i = length(currents1):length(tvec)
        currents1(i) = mu(end);
    end
    movave1 = movmean(currents1, 1000);

    currents2 = [];
    for mu = mus_d
        current = [];
        current(1) = mu;
        for i = 2:(length(tvec) / length(mus_d))
            current(i) = current(i-1) + ((mu - current(i-1)) / tau)*dt + ...
                sigma_d*normrnd(0,1)*sqrt((2*dt) / tau);
        end
        currents2 = horzcat(currents2, current);
    end
    for i = length(currents2):length(tvec)
        currents2(i) = mu(end);
    end
    movave2 = movmean(currents2, 1000);
    
    Iinjs = currents1;
    Iinjd = currents2;
    %Iinjd = zeros(1, length(currents1));
    
    for i = 1:length(tvec)
        if vs(i) == reset
            j = j+1;
            tsum = j*dt;
            vd(i+1) = vd(i) + 10; %elevated by 10 mV
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i)); % 10^-3V * 10^-6 A/V = 10^-9 amps
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vs(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vs(i)));
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT + Iahp(i) + Iinjs(i));  
        elseif vs(i) < theta && j==0
            tsum=0;
            Iahp(i) = 0; %gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vs(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vs(i)));
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));    
       elseif vs(i) < theta && j>0
            j =j+1;
            tsum = j*dt;
            Iahp(i) =gAHP*(EK - vs(i))*exp(-tsum/tauK) ; 
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vs(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vs(i)));
            vs(i+1) = vs(i) + (dt/CS)*((Vrs - vs(i))./RS + (vd(i) - vs(i))./RT +  Iahp(i) + Iinjs(i));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));  
        elseif vs(i) == 10
            j = j+1;
            tsum = j*dt;
            vs(i+1) = reset;
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vs(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vs(i)));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));
        else
            j=1;
            tsum = j*dt;
            vs(i+1) = 10; %for 1ms
            Iahp(i) = gAHP*(EK - vs(i))*exp(-tsum/tauK) ; %activated when fires?
            ICa(i) = gCa*m(i)*h(i)*(Eca -vd(i));
            m(i+1) = m(i) + (dt/taum)*(-m(i) + minf(vs(i)));
            h(i+1) = h(i) + (dt/tauh)*(-h(i) + hinf(vs(i)));
            vd(i+1) = vd(i) + (dt/CD)*((Vrd - vd(i))./RD + (vs(i) - vd(i))./RT + ICa(i) + Iinjd(i));   
        end
    end
    window_size = floor(length(tvec) / length(mus_s)); %step duration of 2s
    num_spikes_s = []; %soma

    cv = [];
    loc_s = cell(size(mus_s));
    isi = cell(size(mus_s));
    for i = 0:length(mus_s)-1
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
ave_freqs_1 = freqs_sum ./ num_iters; %soma 


%plot ICa
figure()
plot(tvec,vs(1:length(tvec)),'k-')

figure()
plot(tvec,ICa)

figure()
plot(tvec,m(1:length(tvec)),'k-')
hold on 
plot(tvec,h(1:length(tvec)),'r-')
hold off

figure()
plot(tvec,vd(1:length(tvec)),'k-')
%% plot figure5
%5A transfer functions
mus_s1= mus_s.*10^3;
figure(1)
plot(mus_s1, ave_freqs_1, 'ok','MarkerFaceColor','k')
hold on
plot(mus_s1, ave_freqs_2, 'ok','MarkerFaceColor','b')
%plot(mus_s1, ave_freqs_3, 'ok','MarkerFaceColor','m')
plot(mus_s1, ave_freqs_4, 'ok','MarkerFaceColor','r')
hold off
ylabel('Mean spike rate (AP/s)')
xlabel('Mean current (pA)')
legend('s','s+d (250 pA)','s+d (750 pA)','Location','Northwest')
title('Modulation of somatic f/I curves')

%5B gain - Done
% slope0 = [];
% for i = 1:length(ave_freqs)
%     g0 = polyfit(pamus(i),ave_freqs(i),1);
%     slope0(i) = g0(1);
% end
% slope0 = max(slope0);

%start from 7th data

slope1 = max(polyfit(mus_s(7:end),ave_freqs_1(7:end),1));
slope2 = max(polyfit(mus_s(7:end),ave_freqs_2(7:end),1));
slope3 = max(polyfit(mus_s(7:end),ave_freqs_3(7:end),1));
slope4 = max(polyfit(mus_s(7:end),ave_freqs_4(7:end),1));
inj_level = [0 0.25 0.5 0.75];
slopes = horzcat(slope1,slope2,slope3,slope4);
slopes_inc = (slopes./slope1)*100; %increase percentage

%fig5Bplotting - done
figure(2)
bar(inj_level,slopes_inc,'k')
xticks([0 0.25 0.5 0.75])
xlabel('Dendritic injection level (nA) n = 10')
ylabel('Gain increase (%)')
title('Gain increase summerized from 10 experiments')

%fig 5C - threshold 

threshold = [];

    for i = 1 : length(ave_freqs_1)
        if ave_freqs_1(i) >= 1
            threshold = horzcat(threshold,mus_s(i));
            break
        else
            i = i+1;
        end
    end
        
     for i = 1 : length(ave_freqs_2)
        if ave_freqs_2(i) >= 1
            threshold = horzcat(threshold,mus_s(i));
            break
        else
            i = i+1;
        end
    end       
    
    for i = 1 : length(ave_freqs_3)
        if ave_freqs_3(i) >= 1
            threshold = horzcat(threshold,mus_s(i));
            break
        else
            i = i+1;
        end
    end

    for i = 1 : length(ave_freqs_4)
        if ave_freqs_4(i) >= 1
            threshold = horzcat(threshold,mus_s(i));
            break
        else
            i = i+1;
        end
    end

figure(3)
plot(inj_level.*1000, threshold.*1000, '-ok','MarkerFaceColor','k')
ylim([ 0 400 ]) 
yticks([0 100 200 300 400])
text(600,380, 'n = 10', 'HorizontalAlignment','center')
ylabel('Mean threshold for somatic current (pA)')
xlabel('Mean dendritic injection currrent (pA)')
title('Mean threshold for somatic current')
    
    
% slope1 = max(polyfit(mus_s(7:end),ave_freqs_1(7:end),1));
% slope2 = max(polyfit(mus_s(5:end),ave_freqs_2(5:end),1));
% slope3 = max(polyfit(mus_s(3:end),ave_freqs_3(3:end),1));
% slope4 = max(polyfit(mus_s(:end),ave_freqs_4(1:end),1));

% slope1 = max(polyfit(mus_s,ave_freqs_1,1));
% slope2 = max(polyfit(mus_s,ave_freqs_2,1));
% slope3 = max(polyfit(mus_s,ave_freqs_3,1));
% slope4 = max(polyfit(mus_s,ave_freqs_4,1));