    clear
    close all
    %% BATCH Analysis ALL videos analysed and stored seperately
    folder = '';
    
    %% New parameter to save the D from all tracks irrespective of the Rsquare
    fitallD = 0;
    % Provide file name base to save output image and .txt data
    ImDataSave = '';
    % Title on output images, specifying treatment or genetic condition
    titleCondition = '';

    % The calibration factor of the microscope, i.e. # pxl = 1 um
    micrCalibration = 1/0.108;

    % The (frame) acquisition time in seconds (s)
    dT = 1;

    % Number of frames to fit for MSD
    frames = 3; % frames-1 = real # of frames

    trackFiles = 'tracks.simple.filtered.txt';
  
    files = dir(folder);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolder = files(dirFlags);
    subFolder(1:2) = [];
    subFolder = {subFolder.name};
    % Read in data from multiple track .txt files belonging to one condition
    numberFiles = length(subFolder);
    for number = 1:numberFiles
        % Compose complete file name
        filenamecomplete = fullfile(folder,subFolder(number),trackFiles);
        filenamecomplete = filenamecomplete{1,1};
        % Read in raw data and put into struct array
        Raw_data(number).data = dlmread(filenamecomplete);
        Raw_data2(number).data = Raw_data(number).data;
        % Conversion of pxl values (x and y) in raw data file(s) into um
        Raw_data(number).data(:,2) = Raw_data(number).data(:,2)/micrCalibration;
        Raw_data(number).data(:,3) = Raw_data(number).data(:,3)/micrCalibration;
    end

    % DATA information on tracks
    % 1st column: time points
    % 2nd column: fitted x-position
    % 3rd column: fitted y-position
    % 4th column: track number true (1) or false (0) label of detected spot
    % 5th column: displacement from the previous time point
    % 6th column: intensity of fitted spot
    % 7th column: fitted sigma of PSF
    % 8th column: fit error

    % Create vector that saves theoretical track length (in # frames for each
    % raw data file) and number of tracks in each file
    trackLength = zeros(1,numberFiles);
    numTracks = zeros(1,numberFiles);

    for file = 1:numberFiles
        sizeRawData = size(Raw_data(file).data);
        trackLength(file) = sizeRawData(1)/(Raw_data(file).data(sizeRawData(1),4) - Raw_data(file).data(1,4) + 1);
        numTracks(file) = Raw_data(file).data(sizeRawData(1),4) - Raw_data(file).data(1,4) + 1;
    end

    FitResults.TotalNumTracks = sum(numTracks);

    % Result for all species together, gives total number of tracks by type (MSD)
    filenamecompleteNumTracks = [folder,ImDataSave,'FitRes_TotalNumTracks.txt'];
    dlmwrite(filenamecompleteNumTracks,FitResults.TotalNumTracks,'delimiter','\t','precision','%.6f','newline','pc')


    for rawDataFile = 1:numberFiles
        % calculate N_TIME_STEPS (for MSD) ~ # frames
        N_TIME_STEPS = zeros(numTracks(rawDataFile),1);
        for trackNumber = 0:numTracks(rawDataFile)-1

            N_TIME_STEPS (trackNumber+1) = sum(Raw_data(rawDataFile).data(:,4) == trackNumber);

        end

        % Making cell array - tracks - for MSD
        tracks = cell(numTracks(rawDataFile), 1);

        for i = 1:numTracks(rawDataFile)
            Time = Raw_data(rawDataFile).data(Raw_data(rawDataFile).data(:,4)== i-1,1); % real time steps    (0 : N_TIME_STEPS(i))' * dT;
            X = Raw_data(rawDataFile).data(Raw_data(rawDataFile).data(:,4)== i-1,2);
            Y = Raw_data(rawDataFile).data(Raw_data(rawDataFile).data(:,4)== i-1,3);
            DetectionTrack = Raw_data(rawDataFile).data(Raw_data(rawDataFile).data(:,4)== i-1,4);
            tracks {i} = [Time X Y DetectionTrack];
        end


        %% MSD direct - calculate for each tau value and diffusion coefficient

        % Number of tau calculated & store MSD values for every tau
        num_tau = max(N_TIME_STEPS)-1;
        nr_track = numTracks(rawDataFile);

        msd = zeros(num_tau,1);

        for dt = 0:num_tau
            if dt == 0
                msd(1,1) = 0;
            end

            sum_of_sq_displ = 0;
            n_sq = 0;
            for Track = 1: nr_track
                % Check if we can get the displacements
                if size(tracks{Track},1) > dt

                    for i = 1:size(tracks{Track},1)-dt
                        % Take into account the gaps
                        if  (tracks{Track}(i+dt,1)-tracks{Track}(i,1)) == dt
                            delta_x = tracks{Track}(i,2)-tracks{Track}(i+dt,2);
                            delta_y = tracks{Track}(i,3)-tracks{Track}(i+dt,3);
                            dd = delta_x^2+delta_y^2;
                            sum_of_sq_displ = sum_of_sq_displ + dd;
                            n_sq = n_sq+1;
                        else
                            for j = 1:(dt-1)
                                if (tracks{Track}(i+dt-j,1)-tracks{Track}(i,1)) == dt
                                    delta_x = tracks{Track}(i,2)-tracks{Track}(i+dt-1,2);
                                    delta_y = tracks{Track}(i,3)-tracks{Track}(i+dt-1,3);
                                    dd = delta_x^2+delta_y^2;
                                    sum_of_sq_displ = sum_of_sq_displ + dd;
                                    n_sq = n_sq+1;
                                    break;
                                end
                            end
                        end
                    end

                end

            end
            msd(dt+1) = (sum_of_sq_displ)/n_sq; % for this particular tau

        end
        % Eliminate Nan values
        msd = msd(~isnan(msd));

        % For the fit/plot of the MSD curve
        time = dT*(0:(length(msd))-1)';

        %     25 % of the average MSD curve
        %     percent = 0.25;
        %     First_25 = 2:(round(percent*length(msd))+1);

        First_25 = 2:size(msd(1:frames),1);
        First_25_msd = msd(First_25);

        % Input for fit function
        xData = time(First_25);
        yData = First_25_msd;

        % Make fit
        [fitresultmsd, gofmsd,output] = createFit(xData, yData);
        coef_MSD_fit = fitresultmsd.p1;

        % Calculate the D = diffusion constant with slope = D*4
        D_msd = coef_MSD_fit/4;

        % fprintf('linear fit of direct MSD curve: D = %d um^2/s \n', D_msd)
        GoodnessFitmsd(rawDataFile) = gofmsd.rsquare;

        figure('Visible','off')
        plot(time,msd,'LineWidth',1)
        hold on
        plot(fitresultmsd, xData, yData)
        title('The direct MSD for different number of tau')
        xlabel('No. of \tau ','Fontsize', 28)
        ylabel(' MSD (\mum^2)','Fontsize', 28)
        set(gca, 'Fontsize', 20);
        data = [time,msd];
        FileName = [folder,subFolder{1,rawDataFile},ImDataSave,'MSD-direct','.txt']
        dlmwrite(FileName,data,'delimiter','\t','precision',3);
%     figure
        %     fitted_values = msd(2:frames)-output.residuals;
        %     scatter(fitted_values,output.residuals)
        %     hline = refline(0);
        %     hline.Color = 'r';

        % File name may also include (raw) data folder
        Filename = [folder,subFolder{1,rawDataFile},ImDataSave,'MSD-direct','.tiff'];

        % Save file (specify image type)
        print('-dtiff',Filename)
        Filename = [folder,subFolder{1,rawDataFile},ImDataSave,'MSD-direct','.eps'];
        print('-depsc','-painters','-loose',Filename)

        %% MSD PER TRACK - Mean MSD
        MSD = zeros(num_tau,nr_track);

        for dt = 0:num_tau
            for Track = 1: nr_track
                if dt ==0
                    MSD(1,Track) = 0;
                else

                    sum_of_sq_displ = 0;
                    n_sq = 0;

                    %       check if we can get displacements for tau=dt
                    if size(tracks{Track},1)>dt

                        for i = 1:size(tracks{Track},1)-dt
                            % Take into account the gaps
                            if  (tracks{Track}(i+dt,1)-tracks{Track}(i,1)) == dt
                                delta_x = tracks{Track}(i,2)-tracks{Track}(i+dt,2);
                                delta_y = tracks{Track}(i,3)-tracks{Track}(i+dt,3);
                                dd = delta_x^2+delta_y^2;
                                sum_of_sq_displ = sum_of_sq_displ + dd;
                                n_sq = n_sq+1;

                            else
                                for j = 1:(dt-1)
                                    if (tracks{Track}(i+dt-j,1)-tracks{Track}(i,1)) == dt
                                        delta_x = tracks{Track}(i,2)-tracks{Track}(i+dt-1,2);
                                        delta_y = tracks{Track}(i,3)-tracks{Track}(i+dt-1,3);
                                        dd = delta_x^2+delta_y^2;
                                        sum_of_sq_displ = sum_of_sq_displ + dd;
                                        n_sq = n_sq+1;
                                        break;
                                    end
                                end
                            end

                        end


                    end
                    MSD(dt+1, Track) = sum_of_sq_displ/n_sq; % for this particular tau
                end
            end
        end

        Sum = nansum(MSD,2);
        NumSum = sum(~isnan(MSD),2);
        AV_MSD = Sum./NumSum;
        AV_MSD = AV_MSD(~isnan(AV_MSD));


        %      25 % of the average MSD curve
        %     percent = 0.25;
        %     FIRST_25 = 2:(round(percent*length(AV_MSD))+1);

        FIRST_25 = 2:size(AV_MSD(1:frames),1);
        FIRST_25_MEAN_MSD = AV_MSD(FIRST_25);
        Time = dT*(0:length(AV_MSD)-1);

        % input for fit function
        xData = Time(FIRST_25)';
        yData = FIRST_25_MEAN_MSD;

        % Make fit
        [fitresultMSD, gofMSD,~] = createFit(xData, yData);
        COEF_MEAN_MSD_FIT = fitresultMSD.p1;

        % slope = D*4, D = diffusion constant
        D_MEAN_MSD_FIT = COEF_MEAN_MSD_FIT/4;

        % fprintf('linear fit of the mean MSD curve: D = %d um^2/s \n', D_MEAN_MSD_FIT)

        GoodnessFitMSD(rawDataFile) = gofMSD.rsquare;
        
        data = [Time',AV_MSD];
        FileName = [folder,subFolder{1,rawDataFile},ImDataSave,'mean_MSD','.txt']
        dlmwrite(FileName,data,'delimiter','\t','precision',3);
        
        figure('Visible','off')
        plot(Time,AV_MSD,'LineWidth',1)
        hold on
        plot(fitresultMSD, xData, yData)
        title('The mean MSD - 2 over different number of tau')
        xlabel('No. of \tau ','Fontsize', 28) % different amount of time steps; 1*tau,2*tau,3*tau
        ylabel('average MSD (\mum^2)','Fontsize', 28)
        set(gca, 'Fontsize', 20);

        %   % plot all MSD curves
        % figure
        % plot(MSD)
        % title('All MSD curves')
        % xlabel('Different # of tau (s)')
        % ylabel('MSD (\mum^2)')

        % File name may also include (raw) data folder
        Filename = [folder,subFolder{1,rawDataFile},ImDataSave,'mean_MSD','.tiff'];

        % Save file (specify image type)
        print('-dtiff',Filename)
        Filename = [folder,subFolder{1,rawDataFile},ImDataSave,'mean_MSD','.eps'];
        print('-depsc','-painters','-loose',Filename)

        %% 3) MSD seperate fit
        % % Store Diffusion coefficients of every track
        ALL_D_COEF = zeros(nr_track,1);
        ALL_R_SQ = zeros(nr_track,1);
        % Make fit on every 25% of the track
        for i = 1:size(MSD,2)
            %             if round(percent*size(MSD(~isnan(MSD(:,i))),1)) > 0
            %             FIRST_25 = 2:(round(percent*size(MSD(~isnan(MSD(:,i))),1))+1);
            if size(MSD(~isnan(MSD(1:frames,i))),1) > 2
                FIRST_25 = 2:size(MSD(~isnan(MSD(1:frames,i))),1);
                TIME = dT*(0:size(MSD(~isnan(MSD(:,i))),1)-1);

                xData = TIME(FIRST_25)';
                yData = MSD(FIRST_25,i);

                % Use function
                [fitresult, gof,~] = createFit(xData, yData);
                % good enough fit
                ALL_R_SQ(i) = gof.rsquare; 
                if gof.rsquare > 0.8||fitallD==1
                    % for every good enough fit the diffusion coefficient calculated by slope = D*4
                    ALL_D_COEF(i) = (fitresult.p1)/4;
                else
                    ALL_D_COEF(i) = 0;
                end
            end
        end

        % only Diff coef of good enough fit
        ALL_D_COEF_filtered = ALL_D_COEF(ALL_D_COEF > 0.0);
        MEAN_D_ALL_TRACKS = mean(ALL_D_COEF_filtered);
        STD_D = std(ALL_D_COEF_filtered);

        % N = # of good enough fit
        N = sum(ALL_D_COEF_filtered~=0);

        %     fprintf('linear fit on every MSD curve: D = %d +- %.3g (mean +- std, N = %d) \n', MEAN_D_ALL_TRACKS , STD_D, N)

        %%

        diffvec = [D_MEAN_MSD_FIT D_msd MEAN_D_ALL_TRACKS STD_D];
        matdiff(rawDataFile,:) = diffvec;

        % 10 log histogram
        ALL_D_LOG = log10(ALL_D_COEF_filtered);
        
         nbins = 80;
         y = ALL_D_LOG;
         [n_sq,xout] = hist(y,nbins);
% 
%         % Data normalisation
         n_sum = sum(n_sq); % number of all data points
         n_norm = n_sq/n_sum; % divide absolute frequency in each bin by total number
% 
        
        % plot
                 figure('Visible','off')


         bar(xout,n_norm); % To plot the histogram
        xlabel('Log diffusion coefficients')
         title('Normalized Bar plot of Log D')
         set(gca, 'Fontsize', 15);

      
 Filename = [folder,subFolder{1,rawDataFile},ImDataSave,'Histogram-D-GaussFit','.tiff'];
% 
%         % Save file (specify image type)
        print('-dtiff',Filename)   
        
        data = [ALL_D_COEF_filtered];
        FileName = [folder,subFolder{1,rawDataFile},ImDataSave,'ALL_D_COEF_filtered','.txt']
        dlmwrite(FileName,data,'delimiter','\t','precision',3);
        %cell2mat(struct2cell(Raw_data(rawDataFile)))
        data = [transpose([0:(nr_track-1)]),ALL_D_COEF,ALL_R_SQ];
        FileName = [folder,subFolder{1,rawDataFile},ImDataSave,'ALL_D_COEF','.txt']
        dlmwrite(FileName,data,'delimiter','\t','precision',3);
        
        %Save tracks.simple with fitted D and Rsquare
        seg_table = array2table((Raw_data2(rawDataFile).data));
      
        track_table = array2table([transpose([0:(nr_track-1)]),ALL_D_COEF,ALL_R_SQ]);
        track_table.Properties.VariableNames = {'Var4';'D';'Rsq'};
      
        seg_table = join(seg_table,track_table,'Keys','Var4');
        data = table2array(seg_table);
        FileName = [folder,subFolder{1,rawDataFile},'\',ImDataSave,'tracks.simple.Dfit','.txt']
        file = fopen(FileName, 'w');
        for ii = 1:size(data, 1)
            fprintf(file, '%0.2f %0.2f %0.2f %0.0f %0.2f %0.2f %0.2f %0.2f  %0.3f  %0.3f\n', data(ii,1),data(ii,2),data(ii,3),data(ii,4),data(ii,5),data(ii,6),data(ii,7),data(ii,8),data(ii,9),data(ii,10));
        end
        dlmwrite(FileName,data,'delimiter','\t','precision',3);
        

    end

      


