function correlation_sum=fista_autocorrelogram(X1_max,clust_index)
si=20;%this is the time per data point
show_ms=100;%prox time window for x axis, 100ms if 100
pad=show_ms*1e3/si;%pad length is data point length
time_res=1;% time resolution for the histogram count
bin=-show_ms:time_res:show_ms;
clust_num=max(clust_index);% Total cluster numbers
map=colormap(jet(clust_num));% Colormap for different clusters
clust=struct();
%% Separate different cluster events out
for i=1:clust_num
    clust(i).LM=X1_max(clust_index==i);
end
p1=zeros(clust_num,clust_num);
p5=zeros(clust_num,clust_num);
p10=zeros(clust_num,clust_num);
correlation_gap=zeros(clust_num,clust_num);
%% Correlogram
% The idea to auto-correlate or cross-correlate is to look at a single event, and count probability distribution of events nearby from
% either the same cluster or other clusters
figure('Unit','Normal','position',[0 0.1 0.6 0.9]);
subplot_num=1;
for j=1:clust_num
    spikes_j=zeros(1,max(X1_max)+2*pad);
    spikes_j(clust(j).LM+pad)=1;
    for i=1:clust_num 
        spikes_i=zeros(1,max(X1_max)+2*pad);
        spikes_i(clust(i).LM+pad)=1;
        dist_prox_index=[];
        for k=1:length(spikes_j)
            %% Look at every spike first
            if spikes_j(k)==1
                % For each spike, use it as a reference, and examine a small spike train within (pad)
                % ms for that spike
                bin_train=spikes_i(k-pad:k+pad);
                % Make the reference spike to 0 so it will not be counted
                bin_train(pad+1)=0;
                % Calculate the distance of other spikes to the reference
                % spike
                dist_prox_index=[dist_prox_index,find(bin_train==1)-pad-1];
            end
        end
        dist_prox_index=dist_prox_index.*si/1e3;
        subplot(clust_num,clust_num,subplot_num)
        if i==j
            map_color=map(j,:);% if auto-correlogram, color-coded by cluster colors
        else
            map_color=[0.3 0.3 0.3];% if cross-correlogram, colored by gray
        end
%         histogram(dist_prox_index,bin,'Normalization','pdf','FaceColor',map_color,'EdgeColor','none')
        histc=histcounts(dist_prox_index,bin);
        histc=histc./length(dist_prox_index);
        bar((bin(1:end-1)+bin(2:end))./2,histc,'FaceColor',map_color,'EdgeColor','none')
        p1(i,j)=round(sum(histc(100:101))./1*show_ms/1,2);
        p5(i,j)=round(sum(histc(96:105))./1*show_ms/5,2);
        p10(i,j)=round(sum(histc(91:110))./1*show_ms/10,2);
        if p1(i,j)<=0.3||(sum([p1(i,j),p5(i,j),p10(i,j)]<=0.5)>1)
            p_color='r';
            correlation_gap(i,j)=1;
        else
            p_color='k';
            correlation_gap(i,j)=0;
        end
        title([num2str(p1(i,j)) '/' num2str(p5(i,j)) '/' num2str(p10(i,j))],'Color',p_color)
        xlim([bin(1),bin(end)])
        y_upper=max([max(histc) 0.02]);
        ylim([0 y_upper])
        subplot_num=subplot_num+1;
    end
end
xlabel('ms')
ylabel('Probability')
correlation_sum=struct('p1',p1,'p5',p5,'p10',p10,'correlation_gap',correlation_gap);
end 
