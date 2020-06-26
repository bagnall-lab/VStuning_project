function process_event_detect(filename_h)
f_mat = dir([filename_h '*.mat']);
for i =1:length(f_mat)
    clearvars S event_index amps poi
     S = load(f_mat(i).name);
     der=struct();
     %% Identify if recordings are for EPSC or EPSP, the "type" descriptions should be put in before event detection
     if isfield(S,'type')
         istype_EPSC=strcmp(S.type{1},'EPSC');
         istype_EPSP=strcmp(S.type{1},'EPSP');
     end
     
    %% If EPSC, a set of parameters are used for detecting EPSC events
    if contains(f_mat(i).name,'EPSC')||istype_EPSC
        % Parameters for EPSCs are defined here
        amp_thre = 6; diff_gap = 240; diff_thre =-8;if_2der=1;event_duration = 640;
        if isfield(S,'Data')
            [der.event_index,der.event_peak,der.amps,der.der_index] = EPSC_detection(S.Data(:,1),S.si,amp_thre,if_2der,diff_gap,diff_thre,event_duration);
        elseif isfield(S,'data_pad')
            [der.event_index,der.event_peak,der.amps,der.der_index] = EPSC_detection(S.data_pad,S.si,amp_thre,if_2der,diff_gap,diff_thre,event_duration);
        else
            warning(['No data in the file: ' f_mat(i).name])
        end
        save(f_mat(i).name,'der','-append');
    end
    %% If EPSP, a different set of parameters are used for detecting EPSC events
    if contains(f_mat(i).name,'EPSP')||istype_EPSP
        % Parameters for EPSPs are defined here
        amp_thre = 0; diff_gap = 480; diff_thre =6;event_duration =1200;
        [der.event_index,der.event_peak,der.amps] = EPSP_detection(S.Data(:,1),S.si,amp_thre,diff_gap,diff_thre,event_duration);
        save(f_mat(i).name,'der','-append');
    end
end
end