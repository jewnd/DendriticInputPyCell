% run for soma + dendrite 250pA data
%% Noisy input current

sigma = 0.3; %nA
%sigma = 0;
tau = 3; %ms

% mus = [0, 0.03, 0.06, 0.09, 0.12, 0.15,...
%     0.18, 0.21, 0.24, 0.27, 0.30, 0.33]; %in nA
% test soma threshold
d_inj = 0.25; %in nA
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

    currents2 = [];
    for mu = mus_d
        current = [];
        current(1) = mu;
        for i = 2:(length(tvec) / length(mus_d))
            current(i) = current(i-1) + ((mu - current(i-1)) / tau)*dt + ...
                sigma*normrnd(0,1)*sqrt((2*dt) / tau);
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
ave_freqs_2 = freqs_sum ./ num_iters; %soma 