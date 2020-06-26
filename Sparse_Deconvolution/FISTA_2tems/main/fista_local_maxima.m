function [X1_max,recon_integral,chemical]=fista_local_maxima(signal,X1,X2,template1,template2,is_plot)
X1_recon=conv(X1,template1);
X2_recon=conv(X2,template2);
peak_width=5;
int_width=50;

X1_max=get_local_maxima_above_threshold(X1,3.5*std(X1),1);
recon_integral=zeros(length(X1_max),2);

% figure;
% subplot(2,1,1)
% hold on;
for i=1:length(X1_max)
    if X1_max(i)+peak_width<length(X1)&&X1_max(i)-peak_width>0
%         plot(X1(X1_max(i)-peak_width:X1_max(i)+peak_width),'k')
        recon_integral(i,1)=sum(X1_recon(X1_max(i)+1:X1_max(i)+int_width))-int_width*max(X1_recon(X1_max(i)+1:X1_max(i)+int_width));
    end
end
% hold off;
% subplot(2,1,2)
% hold on;
for i=1:length(X1_max)
    if X1_max(i)+peak_width<length(X2)&&X1_max(i)-peak_width>0
%         plot(X2(X1_max(i)-peak_width:X1_max(i)+peak_width),'k')
        recon_integral(i,2)=sum(X2_recon(X1_max(i)+1:X1_max(i)+int_width))-int_width*max(X2_recon(X1_max(i)+1:X1_max(i)+int_width));
    end
end
% hold off;

chemical=-recon_integral(:,2)>-recon_integral(:,1);

if is_plot

figure;
plot(signal);
hold on;
scatter(X1_max,signal(X1_max),'r*')
figure;
scatter(-recon_integral(:,1),-recon_integral(:,2),'k');

figure;
subplot(4,1,1)
plot(X1);
hold on;
scatter(X1_max(chemical),X1(X1_max(chemical)),'b')
scatter(X1_max(~chemical),X1(X1_max(~chemical)),'g')
hold off;
subplot(4,1,2)
plot(X2);
hold on;
scatter(X1_max(chemical),X2(X1_max(chemical)),'b')
scatter(X1_max(~chemical),X2(X1_max(~chemical)),'g')
hold off;
subplot(4,1,3)
plot(signal,'k')
hold on;
plot(conv(X1,template1)+conv(X2,template2),'r')
scatter(X1_max(chemical),signal(X1_max(chemical)),'b')
scatter(X1_max(~chemical),signal(X1_max(~chemical)),'g')
hold off;
subplot(4,1,4)
plot(X1_recon,'b')
hold on;
plot(X2_recon,'c')
hold off;
samexaxis('abc','xmt','on','ytac','join','yld',1);
end
end