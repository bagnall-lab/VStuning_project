load('all_traces_padded.mat');
load('EPSC_templates.mat');
Y=smooth(data_pad'-median(data_pad));
% l=6e4:12e4;
%l=1:length(Y);
l=6.48e6:6.6e6;
signal=Y(l,1);
template=fast_EPSC';
opts.backtracking=true;
opts.verbose=false;
opts.lambda=rms(signal).*norm(template).*norminv(0.99);
opts.pos=false;
Xinit=[];
%X = fista_lasso_backtracking_template(Y(l,1), template, Xinit, opts);
%[~,X_alter]=signal_deconv(Y(1:length(Y),1), template,5e4,50,2000);
%% The lasso problem is: argmin( 1/2(Y-A*x1)^2+lambda1*|x1|)
X = fista_lasso_backtracking_template(Y(l,1), template, Xinit, opts);
X_max=get_local_maxima_above_threshold(X,3.5*std(X),0);
%% Plot results
figure;subplot(2,1,1);plot(signal);hold on; scatter(X_max,signal(X_max),'*r');subplot(2,1,2);plot(X);hold on; scatter(X_max,X(X_max),'*r');
samexaxis('abc','xmt','on','ytac','join','yld',1);
