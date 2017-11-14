

clear all;
close all;
numbins=25;
set(0,'defaultAxesFontName', 'Arial');

savepath='/home/ben/Dropbox/5ht3a_project/Dissertation/figures/chapter1/neuron_model/';
filepath='/home/ben/Dropbox/5ht3a_project/neuron_model/MCMC_Circuit_Sept_2017/5ht3optim_PVandPyr_BuildHist_v2.csv';

x=csvread(filepath,1,1);
x=x(:,1:2);
x=x(1:1500,:);
x=x*1000;%convert from uS to nS
for i=length(x(:,1)):-1:1
    if x(i,1)==0&&x(i,2)==0
        x(i,:)=[];
    end
end


[N,C]=hist3(x,[numbins numbins]);
zlim([0 30]);

[val,ind]=max(N);
[maxval,maxind]=max(val);

histcol=floor(find(N==maxval)/length(N(:,1)))+1;
histrow=find(N==maxval)-(histcol-1)*length(N(:,1));
if histrow==0
    histrow=10;
end

disp(strcat('optimal SCpv=',num2str(C{1}(histrow))));
disp(strcat('opimal PVpyr=',num2str(C{2}(histcol))));

figure;
hist3(x,[numbins numbins]);
set(gcf,'renderer','opengl');hold on;
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto')
xlabel('5HT3-Pyr(nS)','FontSize',20);
ylabel('5HT3-PV (nS)','FontSize',20);
title('P(\theta|Pr(Pyr-spike)_e_x_p,Pr(PV-spike)_e_x_p)','FontSize',20);
set(gca,'FontSize',20);
set(gca,'ZTick',0:20:140);
Ztick=cell(length(0:20:140),1);
pp=1;
for i=0:20:140
    Ztick{pp,1}=num2str(round(i/sum(sum(N)),3));
    pp=pp+1;
end
set(gca,'ZTickLabel',Ztick);
saveas(gcf,strcat(savepath,'MCMC_histogram'),'epsc');
saveas(gcf,strcat(savepath,'MCMC_histogram'),'png');


figure;
pcolor(C{1},C{2},N');hold on;

% for i=1:length(x(:,1))
%     plot(x(1:i,1),x(1:i,2),'-om','LineWidth',2);hold on;
%     pause(.1);
%     
% end

set(gca,'FontSize',25,'Linewidth',3);
xlabel('5HT3-Pyr (nS)','FontSize',25);
ylabel('5HT3-PV (nS)','FontSize',25);
axis square;
figure;%set(gcf,'Renderer','painters');
pcolor(C{1},C{2},N');hold on;
% F(150) = struct('cdata',[],'colormap',[]);
% vidObj = VideoWriter('/home/ben/Dropbox/MCMC_randwalk.avi');
% open(vidObj);
 plot(x(1:end,1)+.0005,x(1:end,2)+.001,'-m','LineWidth',3);hold on;
 plot(x(1,1),x(1,2),'ok','MarkerFaceColor','k');hold on;

set(gca,'FontSize',25,'Linewidth',3);hold on;
xlabel('5HT3-Pyr (\muS)','FontSize',25);hold on;
ylabel('5HT3-PV (\muS)','FontSize',25);hold on;
axis square;
% currFrame=getframe;
% writeVideo(vidObj,currFrame);



% xlabel('5HT3-Pyr (\muS)','FontSize',25);hold on;
% ylabel('5HT3-PV (\muS)','FontSize',25);hold on;
xlabel('5HT3-Pyr (nS)','FontSize',25);hold on;
ylabel('5HT3-PV (nS)','FontSize',25);hold on;
axis square;
saveas(gcf,strcat(savepath,'MCMC_randomwalk'),'epsc');
saveas(gcf,strcat(savepath,'MCMC_randomwalk'),'png');



pvpyrratio=x(:,2)./x(:,1);
figure;hold on;
hist(pvpyrratio,20);hold on;
% xlim([0 10]);
set(gca,'LineWidth',2,'FontSize',20);
xlabel('5HT3-PV/5HT3-Pyr','FontSize',20);
ylabel('Counts','FontSize',20);
box off;
axis square;
saveas(gcf,strcat(savepath,'PV_Pyr_Ratio_Histogram'),'epsc');
saveas(gcf,strcat(savepath,'PV_Pyr_Ratio_Histogram'),'png');
% xax = 0 : .005 : .06; % axis x, which you want to see
% yax = 0 : .05 :1; % axis y, which you want to see
% 
% [X,Y] = meshgrid(xax,yax); % important for "surf" - makes defined grid
% 
% pdf = hist3([x(:,1) , x(:,2)],{xax yax}); % standard hist3 (calculated for yours axis)
% pdf_normalize = (pdf' ./ length(x(:,1))); % normalization means devide it by length of
% % data_x (or data_y)
% figure()
% surf(X,Y,pdf_normalize)


%regress values outputted from MCMC
% figure;plot(x(1,1),x(1,2),'ok','MarkerFaceColor','k');hold on;xlim([0 .12]);ylim([0 1.8]);
% for i=2:length(x(:,1))
%     hold on;plot([x(i-1,1) x(i,1)],[x(i-1,2) x(i,2)],'k');pause(0.1);hold on;
% end
