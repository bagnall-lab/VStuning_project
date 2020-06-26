m=5;l=100;
template=rand(m,1);
E=zeros(l,1);
E_alt=E;
% for i=0:l
% %signal=rand(n,1);
% n=m+i;
% t_m=toeplitz([template;zeros(n-m+1,1)],[template(1);zeros(n-m+1,1)]);
% coef=t_m'*t_m;
% if i==0
%     E_0=max(eig(coef));
% else
%     E(i)=max(eig(coef));
% end
% end

n=m+l;
t_m=toeplitz([template;zeros(n-m+1,1)],[template(1);zeros(n-m+1,1)]);
coef=t_m'*t_m;
figure;heatmap(coef)
for i=0:l-1
    E_alt(l-i)=max(eig(coef(1:end-i,1:end-i)));
end
figure;plot(diff(E_alt))
%E_alt(end)./E_alt(5)
%f=fit([1:length(E)]',E./E(1),'exp1')