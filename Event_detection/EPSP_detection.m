function [event_index,event_peak,amps] = EPSP_detection(W,si,amp_thre,diff_gap,diff_thre,event_duration)
%% calculate the difference with 240us as "1st derivative" to detect event
diff_gap = diff_gap/si;
event_duration =event_duration/si;
data_s = smooth(W);  %smooth the data
diff_ = data_s(1+diff_gap:end)-data_s(1:end-diff_gap)-diff_thre; % caculate the difference between 240us
crossing_ = diff_(1:end-1).*diff_(2:end)<0;% find the crossing for -8pA
r_index = find(crossing_.*(diff_(1:end-1)>0)); % find the on-rise crossing point
f_index = find(crossing_.*(diff_(1:end-1)<0)); % find the on-fall crossing point
if length(r_index) ~= length(f_index)
    warning('rising and falling not match')
    if length(r_index) == length(f_index)+1
        r_index(end)=[];
    elseif length(f_index) == length(r_index)+1
        f_index(1)=[];
    end
end

%% calculate amplitude for each events detected
raw_index=r_index;
raw_l = length(raw_index);
amp_raw = zeros(raw_l,1);
peak_index = zeros(raw_l,1);
for i = 1:raw_l
    if i == raw_l
        duration = min(event_duration,length(data_s)-raw_index(i)); % when it is the last EPSC, make sure duration of EPSC didn't go beyond the end of the whole trace
    else
        duration = min(raw_index(i+1)-raw_index(i),32);% Duration of EPSC is the smaller one of either 640us or before next EPSC comes
    end
    [peak_value,peak_index(i)] = max(data_s(raw_index(i):raw_index(i)+duration));
    amp_raw(i)=peak_value-data_s(raw_index(i));
    peak_index(i)=raw_index(i)+peak_index(i)-1;
    % Calculate the amplitude of each EPSC, here the algorithm is:
    % max value within event_duration - value at the index point
end

%% set amplitudes threshold for EPSC detection
amps = amp_raw(amp_raw>amp_thre);
event_index = r_index(amp_raw>amp_thre);
% der_index is the second column of all_index, which keeps track of the event types (1der or 2der) 
event_peak=peak_index(amp_raw>amp_thre);
end