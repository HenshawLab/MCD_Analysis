close all
clc


%% Main 

if RUN_PARTICLELOCATION == true
    for iB = 1:NBio
        
        biostr = ['Bio' num2str(BioReps(iB))];
    
        for iR = 1:NRep
    
            bacteria = cell(NLoop,NPos);
            repstr = ['Rep' sprintf('%02d',Reps(iR))];
    
            % Load preprocessing
            load([Output_Background ExpName '_' biostr '_' repstr '_MeanBackground.mat']);
            load([Output_Pretrack ExpName '_' biostr '_' repstr '_PretrackingParams.mat']);
            load([Output_Cropping ExpName '_' biostr '_' repstr '_CroppingLimits.mat']);
                        
            for ip = 1:NPos
                
                concstr = ['C' num2str(ip-1)];

                % Get image list
                ImgDir = [MainDir biostr '/' repstr '/' concstr '/'];
                if strcmp(filenaming,'custom')
                    for i = 1:NLoop
                        imglist(i).name = ['Position ' num2str(ip) '_t' sprintf('%02d',i-1) '_ch00' imgextension];
                    end
                else
                    cd(ImgDir); imglist = dir(imgextension); cd(WorkingDir);
                end

                for iL = 1:NLoop

                    img = imread([ImgDir imglist(iL).name]);              
                    if size(img,3) > 1
                        img = double(rgb2gray(img));
                    else
                        img = double(img);
                    end
                    
                    % Invert (if appropriate) and remove background
                    img = mat2gray(img) - Bimg(:,:,ip);
                    if strcmp(particle_type,'dark')
                        img = imcomplement(img);
                    end                    

                    % Crop image
                    img = img(yrange(ip,1):yrange(ip,2),:);

                    % Find particles, convert to microns
                    b = bpass(img.*255,BPASS(ip,1),BPASS(ip,2));
                    pk = pkfnd(b,PKFND(ip,1),PKFND(ip,2));
                    cnt = cntrd(b,pk,CNT(ip));

                    % Store data
                    bacteria{iL,ip} = cnt(:,[1,2]).*PixToMum;
                    clear b pk cnt

                end % End of looping over frames
    
            end % End of looping over positions (NPos)

                params.repstr = repstr;
                params.biostr = biostr;
                params.BPASS = BPASS;
                params.PKFND = PKFND;
                params.CNT = CNT;
                params.particle_type = particle_type;
                params.PixToMum = PixToMum;
                params.NLoop = NLoop;
                params.BioReps = BioReps;
                params.Reps = Reps;
                params.ExpName = ExpName;
                params.NPos = NPos;
                params.chem_location = chem_location;

                save([ExperimentOutDir 'Bacteria_Positions/' biostr '_' repstr '_PositionData.mat'],...
                    'bacteria','params');

                clear bacteria params
            
        end % End of technical replicates (NRep)
        
    end % End of biological replicates (NBio)

end % End of pretracking loop

%% Main - analysis

if RUN_ANALYSIS == true
    
    % Preallocate
    Nupper = zeros(NLoop,NPos,NRep,NBio); Nlower = Nupper;
    for ip = 1:NPos
        beta_output(ip).beta = [];
        beta_output(ip).time = [];
        beta_output(ip).BioRep = [];
        beta_output(ip).TechRep = [];
    end
    beta_all = []; time_all = [];

    for iB = 1:NBio

        biostr = ['Bio' num2str(BioReps(iB))];

        for iR = 1:NRep

            repstr = ['Rep' sprintf('%02d',Reps(iR))];
            % Load pretracking and cropping limits
            load([ExperimentOutDir 'Bacteria_Positions/' biostr '_' repstr '_PositionData.mat'])
            load([Output_Cropping ExpName '_' biostr '_' repstr '_CroppingLimits.mat']);

            t0 = [0:NLoop-1].*imaging_period;
        
            for ip = 1:NPos

                Channel_Width = abs((yrange(ip,2)-yrange(ip,1))).*PixToMum;
                NBins = ceil(Channel_Width/BinW);
                BinEdges = [0:NBins].*BinW;

                for ic = 1:length(BinEdges)-1
                    centres(ic) = (BinEdges(ic+1)+BinEdges(ic))./2;
                end

                rhomap = zeros(length(BinEdges)-1,NLoop,NPos);

                % Counting and heatmaps
                for it = 1:NLoop
                    
                    % Get y positions for this time
                    y = bacteria{it,ip}; y = y(:,2);

                    % Count particles in top and bottom
                    ytop = y(y<=accum_width);
                    ybottom = y(y>=(Channel_Width-accum_width));

                    nu(it) = numel(ytop); nl(it) = numel(ybottom);

                    h = histogram(y,BinEdges,'normalization','pdf');
                    rhomap(:,it,ip) = h.Values';
%                     
%                     % For saving purposes
                    Nupper(it,ip,iR,iB) = nu(it);
                    Nlower(it,ip,iR,iB) = nl(it);

                end % End of looping over time

                % Beta calculation
                beta = (nu - nl)./(nu(end)+nl(end));
                if strcmp(chem_location,'bottom')
                    beta = -1.*beta;
                end

                time = t0 + (ip-1)*imaging_period;

                beta_output(ip).beta = [beta_output(ip).beta;beta];
                beta_output(ip).time = [beta_output(ip).time;time];
                beta_output(ip).BioRep = [beta_output(ip).BioRep,iB];
                beta_output(ip).TechRep = [beta_output(ip).TechRep,iR];

                beta_all = [beta_all,beta'];
                time_all = [time_all,time'];         

                close all
                % Heatmap
                fh = figure;
                set(fh,'color','white');
                imagesc(time./60,centres,squeeze(rhomap(:,:,ip)));
                colormap('hot');
                set(gca,'YDir','reverse');
                ylabel('Position y \mum','Interpreter','LaTex');
                xlabel('Time t min','Interpreter','LaTex');
                set(gca,'fontsize',10);
                set(gca,'XLim',[min(time./60),max(time./60)],'YLim',[min(centres),max(centres)]);
                title([biostr ' ' repstr ' ' 'Channel ' num2str(ip)]);
                colorbar;
                saveas(gcf,[FigDir 'Heatmaps/' ExpName '_' biostr '_' repstr '_Channel-' num2str(ip) '.fig'])
                saveas(gcf,[PNGDir 'Heatmaps/' ExpName '_' biostr '_' repstr '_Channel-' num2str(ip) '.png'])
                close(gcf);

            end % End of looping over positions

        end % End of technical replicates (NRep)

    end % End of biological reps (NBio)

    save([ExperimentOutDir ExpName '_beta_analysis.mat'],...
        'accum_width','beta_all','time_all','Nupper','Nlower',...
        'beta_output');

    save([BetaDir ExpName '_beta_analysis.mat'],...
        'accum_width','beta_all','time_all','Nupper','Nlower',...
        'beta_output');

end % End of analysis loop

