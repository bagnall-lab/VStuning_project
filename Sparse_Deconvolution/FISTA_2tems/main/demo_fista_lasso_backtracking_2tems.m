function process_fista_lasso(filename)
S=load(filename);
T=load('EPSC_templates.mat');
Y=smooth(S.data_pad'-median(S.data_pad));
% l=6e4:12e4;
l=1:length(Y);
%l=6.48e6:6.6e6;
signal=Y(l,1);
EPSC_w1=T.fast_EPSC(1:441)';
EPSC_w2=T.slow_EPSC';
template1=EPSC_w1;
%template2=fast_EPSC(1:length(slow_EPSC))'-slow_EPSC';
alpha=EPSC_w1'*EPSC_w1/(EPSC_w1'*EPSC_w2);
template2=EPSC_w1-alpha.*EPSC_w2;
opts.backtracking=true;
opts.verbose=true;
opts.lambda1=rms(signal).*norm(template1);
opts.lambda2=rms(signal).*norm(template2);
opts.pos=false;
Xinit=[];

%[~,Xinit]=signal_deconv(Y(1:length(Y)/100,1), template,5e4,50,2000);
%X = fista_lasso_backtracking_template(Y(l,1), template, Xinit, opts);
[X1,X2,cost_matrix] = fista_lasso_backtracking_2tems(Y(l,1), template1,template2, Xinit,Xinit, opts);

%% The lasso problem is: argmin( 1/2(Y-A*x1-B*x2)^2+lambda1*|x1|+lambda2*|x2|)

end
