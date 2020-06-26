function process_fista_main(filename)

%% Events detection using a two-template sparse deconvolution method(fista_2tem)
process_fista_2tem_detect(filename);
%process_fista_lasso(filename)
%% Clustering events
process_fista_clust(filename);
%% Calculate the amplitude of events 
process_fista_get_amps(filename)
%% Plot the results from fista and clustering
process_plot_fista_results(filename);

end
