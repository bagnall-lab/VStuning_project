%{
  This funtion is used to convert .abf to .mat
  The .mat file stores all data in a structure with fields like:

  name: name of the original file
  Data: recording data
  si: Sampling intervals
  header: header of the file
%}

function process_abf2mat(filename_h)
f_abf = dir([filename_h '*.abf']);
for i =1:length(f_abf)
    clearvars name Data si header
    [~,name,~] = fileparts(f_abf(i).name);
    [Data,si,header] = abfload(f_abf(i).name);
    switch type_of_rec(header,Data(:,:,1))
        case 1
            type={'EPSC'};
        case 2
            type={'IPSC'};
        case 3
            type={'IPSP'};
        case 4
            type={'EPSP'};
    end
    save([name '.mat'],'name','Data','si','header','type');
end
end

function rec_type = type_of_rec(h,data)
    rec_type = 0;
    if (strcmp(h.recChUnits{1},'pA'))
        if median(data(:,1))>mean(data(:,1)) % spikes are down, EPSC
            rec_type = 1;
        else % spikes are up, IPSC
            rec_type = 2;
        end
    elseif (strcmp(h.recChUnits{1},'mV'))
        if median(data(:,1))>mean(data(:,1)) % spikes are down, IPSP
            rec_type = 3;
        else % spikes are up, EPSP
            rec_type = 4;
        end 
    end
end