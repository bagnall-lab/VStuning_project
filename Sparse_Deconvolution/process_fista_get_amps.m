function process_fista_get_amps(filename)
    S=load(filename);
    fista=S.fista;
    [fista.amps,fista.amp_index]=get_amplitude(S.data_pad',S.fista.X1_max);
    save(filename,'fista','-append')
end

function [amp_raw,peak_index]=get_amplitude(data_s,raw_index)
raw_l = length(raw_index);
amp_raw = zeros(raw_l,1);
peak_index = zeros(raw_l,1);
event_duration=32;
for i = 1:raw_l
    if i == raw_l
        duration = min(event_duration,length(data_s)-raw_index(i)); % when it is the last EPSC, make sure duration of EPSC didn't go beyond the end of the whole trace
    else
        duration = min(raw_index(i+1)-raw_index(i),event_duration);% Duration of EPSC is the smaller one of either 640us or before next EPSC comes
    end
    [peak_value,peak_index(i)] = min(data_s(raw_index(i):raw_index(i)+duration));
    amp_raw(i)=data_s(raw_index(i))-peak_value;
    peak_index(i)=raw_index(i)+peak_index(i)-1;
    % Calculate the amplitude of each EPSC:
    % max value within event_duration - value at the index point
end
end