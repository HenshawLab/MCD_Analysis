%% Plot each biological replicate separately
close all

load([ExperimentOutDir ExpName '_beta_analysis.mat']);
repcol = lines(6);

for iB = 1:NBio

    fh = figure;
    set(fh,'color','white');
    box on; hold on;
    
    for ip = 1:NPos

        subplot(2,3,ip)
        t = beta_output(ip).time;
        beta = beta_output(ip).beta;
    
        for iR = 1:NRep
            ph(iR) = plot(t(iR,:)./60,beta(iR,:),'col',repcol(iR,:),'linewidth',1.5); hold on;
            names{iR} = ['Rep ' num2str(Reps(iR))];
        end

        
    end % End of looping over channels

    for ip = 1:NPos
        subplot(2,3,ip);
        xlabel('Time t min'); ylabel('\beta(t)');
        title(['C = ' num2str(concentrations(ip)) '\muM']);
        set(gca,'fontsize',12);
        ax = gca; ax.LineWidth = 1;
        xl = get(gca,'XLim');
        plot(xl,[0,0],'--k','linewidth',0.5,'HandleVisibility','off');
        xlim(xl);

        if ip == 1
            leg = legend(ph,names,'location','northeast');
        end
    end
    sgtitle(ExpName);
    set(gcf,'Position',get(0,'ScreenSize'));

    saveas(gcf,[FigDir ExpName '_BetaCurves_Bio' num2str(BioReps(iB)) '.fig']);
    saveas(gcf,[PNGDir ExpName '_BetaCurves_Bio' num2str(BioReps(iB)) '.png']);

    saveas(gcf,[BetaFigDir ExpName '_BetaCurves_Bio' num2str(BioReps(iB)) '.fig']);
    saveas(gcf,[BetaPNGDir ExpName '_BetaCurves_Bio' num2str(BioReps(iB)) '.png']);

end

close all

%% Plot everything on one figure with SE shading
% 
fh = figure;
set(fh,'color','white'); box on; hold on;
betacols = plasma(NPos);

for ip = 1:NPos

    t = (beta_output(ip).time)./60;
    beta = beta_output(ip).beta;

    pP(ip) = plot(t(1,:),mean(beta,1),'linewidth',1.5,'color',betacols(ip,:)); hold on;

    ppatch(ip) = errpatch(t(1,:),mean(beta,1),std(beta,[],1).*size(beta,1),betacols(ip,:),0.3);

    cnames{ip} = ['C = ' num2str(concentrations(ip)) '\muM'];

end

xlabel('Time t min'); ylabel('\beta(t)');
xl = get(gca,'XLim');
plot(xl,[0,0],'--k','linewidth',0.5,'HandleVisibility','off');
xlim(xl);
ax = gca; ax.LineWidth = 1;
set(gcf,'Position',get(0,'ScreenSize'));
set(gca,'fontsize',15);
leg = legend(pP,cnames,'location','southeast');
set(leg,'box','off');
axis square

title(ExpName)

saveas(gcf,[FigDir ExpName '_BetaCurves.fig']);
saveas(gcf,[PNGDir ExpName '_BetaCurves.png']);

saveas(gcf,[BetaFigDir ExpName '_BetaCurves.fig']);
saveas(gcf,[BetaPNGDir ExpName '_BetaCurves.png']);

