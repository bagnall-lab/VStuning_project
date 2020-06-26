function process_fista_clust(filename)
% Cluster the events with isosplit5 algorithm, also make sure the labels 
% will follow them after events are separated into different files.
    S=load(filename);
    fista=S.fista;
    data_s=smooth(S.data_pad);
    %% Setting up several parameters for clustering
    %X1_max: onset time of each fast event
    
    %X1_prox: deconvoluted signals around X1_max, -5:5
    
    %X1_waveforms: original signal of each event, from onset (0) to 70 (1.4ms)
    
    %X1_integral: the sum value of deconvoluted signals during X1_prox,
    %can be used to estimate the amplitude of EPSC
    
    fista.X1_max=fista.X1_max(fista.X1_max>5&fista.X1_max<(length(fista.X1)-5));
    fista.X1_prox=zeros(length(fista.X1_max),11);
    fista.X1_waveform=zeros(length(fista.X1_max),71);
    fista.X1_integral=zeros(length(fista.X1_max),1);
    fista.X1_std=zeros(length(fista.X1_max),1);
    fista.X12_ratio=sum(abs(fista.X1))./sum(abs(fista.X2));
    
    %% Each parameter is defined in the loop here 
    for i=1:length(fista.X1_max)
        fista.X1_prox(i,:)=fista.X1(fista.X1_max(i)-5:fista.X1_max(i)+5)';
        fista.X1_waveform(i,:)=data_s(fista.X1_max(i):fista.X1_max(i)+70)'-data_s(fista.X1_max(i));
        fista.X1_integral(i)=sum(fista.X1(fista.X1_max(i)-5:fista.X1_max(i)+5));
        fista.X1_std(i)=std(fista.X1(fista.X1_max(i)-5:fista.X1_max(i)+5)); 
        %fista.X1_neighbor(i,:)=fista.X1(fista.X1_max(i)-5:fista.X1_max(i)+50)';
    end
    
    opts.isocut_threshold=1;% setting threshold for clustering, see isosplit5 for specifics
    %% Main clsutering algorithm, in theory, either X1_prox, X1_waveform, X1_integral can be used for clustering, depending on what aspect of EPSCs is most distinctive.
    % In practice, we found X1_prox usually performs better
    fista.X1_clust=isosplit5(fista.X1_prox',opts)';
    %% clustering result also saved in fista -> X1_clust
    save(filename,'fista','-append')
end