%% Sample starting criteria (assuming 6 channels)
BPASS = [1,7;... % 1
    1,7;... % 2
    1,7;... % 3
    1,7;... % 4
    1,7;... % 5
    1,7]; % 6
PKFND = [5,7;... % 1
    5,7;... % 2
    5,7;... % 3
    5,7;... % 4
    5,7;... % 5
    5,7]; % 6

CNT = [7 7 7 7 7 7];

%% Make background images

if RUN_BACKGROUND == true
    for iB = 1:NBio
    
        biostr = ['Bio' num2str(BioReps(iB))];
    
        for iR = 1:NRep
    
            repstr = ['Rep' sprintf('%02d',Reps(iR))];
            
            % Loop over positions, produce a NxMxNLoop matrix for the
            % background
            for ip = 1:NPos
    
                concstr = ['C' num2str(ip-1)];
    
                % Get list of images, setup background preallocation
                ImgDir = [MainDir biostr '/' repstr '/' concstr '/'];
                if strcmp(filenaming,'custom')
                    for i = 1:NLoop
                        imglist(i).name = ['Position ' num2str(ip) '_t' sprintf('%02d',i) '_ch00' imgextension];
                    end
                else
                    cd(ImgDir); imglist = dir(imgextension); cd(WorkingDir);
                end
                img = imread([ImgDir imglist(1).name]);
                
                % Preallocate background
                if ip == 1
                    Bimg = zeros(size(img,1),size(img,2),NPos);
                end
    
                for iL = 1:NLoop
                    
                    % Check if a 3 channel (RGB) image, correct if so
                    img = imread([ImgDir imglist(iL).name]);
                    if size(img,3) > 1
                        img = double(rgb2gray(img));
                    else
                        img = double(img);
                    end
    
                    Bimg(:,:,ip) = Bimg(:,:,ip) + mat2gray(img);
    
                end % End of looping over frames (NLoop)
    
            end % End of looping over positions (NPos)
    
            Bimg = Bimg./NLoop;
            save([Output_Background ExpName '_' biostr '_' repstr '_MeanBackground.mat'],'Bimg');
    
        end % End of looping over technical reps (NRep)
    
    end % End of looping over biologicla reps (NBio)

end % End of if statement

%% Cropping

if RUN_CROPPING == true
    for iB = 1:NBio
    
        biostr = ['Bio' num2str(BioReps(iB))];
    
        for iR = 1:NRep
    
            repstr = ['Rep' sprintf('%02d',Reps(iR))];
            
            % Loop over positions, user defined inputs for cropping limits
            figure; hold on; set(gcf,'Position',get(0,'ScreenSize'));
    
            xrange = zeros(NPos,2); yrange = xrange;
    
            for ip = 1:NPos
    
                concstr = ['C' num2str(ip-1)];
    
                % Get list of images
                ImgDir = [MainDir biostr '/' repstr '/' concstr '/'];
                if strcmp(filenaming,'custom')
                    for i = 1:NLoop
                        imglist(i).name = ['Position ' num2str(ip) '_t' sprintf('%02d',i) '_ch00' imgextension];
                    end
                else
                    cd(ImgDir); imglist = dir(imgextension); cd(WorkingDir);
                end
    
                % Read image, correct to grayscale
                img = imread([ImgDir imglist(1).name]);
                if size(img,3) > 1
                    img = double(rgb2gray(img));
                else
                    img = double(img);
                end
    
                % Plot and ask for limits.
                imshow(img,[]); hold on; 
                title([ExpName ' ' biostr ' ' repstr ' ' concstr],...
                    'Interpreter','None');
    
                disp('Two clicks, top then bottom, at channel boundaries');
                [xi,yi] = ginput(2);
    
                xrange(ip,:) = round(xi');
                yrange(ip,:) = round(yi');
    
            end % End of looping over positions (NPos)
            close all
    
            save([Output_Cropping ExpName '_' biostr '_' repstr '_CroppingLimits.mat'],...
                'xrange','yrange');
            clear xrange yrange
    
        end % End of looping over technical reps (NRep)
    
    end % End of looping over biological reps (NBio)
end % End of if statement
%% Pretracking parameters
close all
if RUN_PRETRACKPARAMETERS == true
    for iB = 1:NBio
    
        biostr = ['Bio' num2str(BioReps(iB))];
    
        for iR = 1:NRep
    
            repstr = ['Rep' sprintf('%02d',Reps(iR))];
            
            % Load background and cropping limits
            load([Output_Background ExpName '_' biostr '_' repstr '_MeanBackground.mat']);
            load([Output_Cropping ExpName '_' biostr '_' repstr '_CroppingLimits.mat'])
            
            figure; hold on; set(gcf,'Position',get(0,'ScreenSize'));

            for ip = 1:NPos
    
                concstr = ['C' num2str(ip-1)];
    
                % Get list of images, setup background preallocation
                ImgDir = [MainDir biostr '/' repstr '/' concstr '/'];
                if strcmp(filenaming,'custom')
                    for i = 1:NLoop
                        imglist(i).name = ['Position ' num2str(ip) '_t' sprintf('%02d',i-1) '_ch00' imgextension];
                    end
                else
                    cd(ImgDir); imglist = dir(imgextension); cd(WorkingDir);
                end

                % Take a mid-time image where particles should be dispersed
                img = imread([ImgDir imglist(floor(NLoop/2)).name]);              
                if size(img,3) > 1
                    img = double(rgb2gray(img));
                else
                    img = double(img);
                end

                % Remove background and plot
                img = mat2gray(img) - Bimg(:,:,ip);
                % Invert if particles are dark
                if strcmp(particle_type,'dark')
                    img = imcomplement(img);
                end
                imshow(img,[]); hold on;

                HAPPY = false;

                while HAPPY ~= true                  
                    
                    % Find particles
                    b = bpass(img.*255,BPASS(ip,1),BPASS(ip,2));
                    pk = pkfnd(b,PKFND(ip,1),PKFND(ip,2));
                    cnt = cntrd(b,pk,CNT(ip));

                    pts = plot(cnt(:,1),cnt(:,2),'rx','linewidth',1);

                    HAPPY = input('Happy? 0 = No, 1 = Yes');

                    if HAPPY == 0

                        disp(['BPASS: ' num2str(BPASS(ip,1)) ' ' num2str(BPASS(ip,2))]);
                        BPASS(ip,:) = input('Input new BPASS: [bp1,bp2] =   ');
                        disp(['PKFND: ' num2str(PKFND(ip,1)) ' ' num2str(PKFND(ip,2))]);
                        PKFND(ip,:) = input('Input new PKFND: [pk1 pk2] =  ');
                        disp(['CNT:' num2str(CNT(ip))]);
                        CNT(ip) = input('New CNT: CNT =  ');
                        delete(pts); clc;

                    else
                        HAPPY = 1;
                        break
                    end % End of checking 

                end % End of while loop
    
            end % End of looping over positions (NPos)
            close all
    
            save([Output_Pretrack ExpName '_' biostr '_' repstr '_PretrackingParams.mat'],...
                'BPASS','PKFND','CNT');
%             clear BPASS PKFND CNT
    
        end % End of looping over technical reps (NRep)
    
    end % End of looping over biological reps (NBio)
end % End of if statement