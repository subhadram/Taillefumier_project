ampp = 0.5:0.5:2.5
%ampp = [1.5];
psi = 0:0.1:2.0 ;
psi = pi*psi;
%freq = [6,7,8,9,10];
freq = [10];
freq = 2*pi/1000.0*freq;
w00 = [0.2,0.4,0.8,1.6,2.5];
%w00 = [2.5];
taut = [5,10,20,30,40]
load('ang10.mat')
load('angsem10.mat')

for t=1:length(taut)
    figure(t)
    %title(sprintf('tau = %d msec 10 Hz',taut(t)))
    for w = 1:length(w00)
        errorbar(ampp,squeeze(ang10(t,w,:)),squeeze(angsem10(t,w,:)),'-o')
        hold on;
    end
    ylim([0,2*pi])
    title(sprintf('tau = %d msec 10 Hz',taut(t)))
    legend('recurrent weight = 0.2','recurrent weight = 0.4','recurrent weight = 0.8','recurrent weight = 1.6','recurrent weight = 2.5')
    xlabel('Coupling strength')
    ylabel('Optimal phase difference')
   set(gca,'fontsize',14)
end

for t=1:length(taut)
    figure(t+5)
    %[X,Y]=meshgrid(w00,ampp)
     %surf(w00,ampp,squeeze(ang6(t,:,:)))
     %view(2)
     colormap(jet)
     imagesc(squeeze(ang10(t,:,:)))
     xticks(1:1:5)
     yticks(1:1:5)
        %errorbar(ampp,squeeze(ang6(t,w,:)),squeeze(ang6sem(t,w,:)),'-o')
    yticklabels({'0.2','0.4','0.8','1.6','2.5'})
    xticklabels({'0.5','1.0','1.5','2.0','2.5'})
    colorbar
    title(sprintf('tau = %d msec 10 Hz',taut(t)))
    %legend('recurrent weight = 0.2','recurrent weight = 0.4','recurrent weight = 0.8','recurrent weight = 1.6','recurrent weight = 2.5')
    ylabel('Recurrent strength')
    xlabel('Coupling strength')
    h=colorbar
    caxis([0 6.3])
    %set(gca,'YTick',ampp)
    ylabel(h, 'Optimal phase difference')
    set(gca,'fontsize',14)
end
