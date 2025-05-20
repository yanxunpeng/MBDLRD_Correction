function MBDLRD_Correction()
    clc;    close all;
%------------load data
    load 'seafloorWithErrors.mat' raster_error raster_R;
    [row,col] = size(raster_error);
    cellWid=raster_R.CellExtentInWorldX;%grid width
    raster_correct = raster_error;
    % mesh(raster_error);

    load 'realSeafloor.mat' raster_real;
    % mesh(raster_real);
%------------Swath identification, segmentation, and reconstruction
    swathCount=4;
    Xcoor=(cellWid:cellWid:cellWid*row)';
    surfaceSoundSpeed=1500;
    minDepthRaster = min(raster_error, [], 'all', 'omitnan') ;
    swath=cell(col,swathCount);
    load SwathFlag.mat SwathFlag;
    for i=1:1:col
        % disp(i);
        terrainProfile=raster_error(:,i);
        for j=1:swathCount%1-5,2-9,6-13,10-14
            if j==1
                startIndex=1;
                endIndex=5;
            elseif j==2
                startIndex=2;
                endIndex=9;
            elseif j==3
                startIndex=6;
                endIndex=13;
            elseif j==4
                startIndex=10;
                endIndex=14;
            end
            SwathStartIndex=SwathFlag(i,startIndex);
            SwathEndIndex=SwathFlag(i,endIndex);
            swath{i,j}.startIndex = SwathStartIndex;
            swath{i,j}.endIndex = SwathEndIndex;
            swath{i,j}.depth=terrainProfile(SwathStartIndex:SwathEndIndex);
            swath{i,j}.x = Xcoor(SwathStartIndex:SwathEndIndex);
            swath{i,j}.across0=swath{i,j}.x-swath{i,j}.x(1);%0开始
            swath{i,j}.across_=swath{i,j}.across0-swath{i,j}.across0(end)/2;%-负号开始
            swath{i,j}.travelTime=sqrt((swath{i,j}.depth).^2+(swath{i,j}.across_).^2)/surfaceSoundSpeed;%单程旅行时
            swath{i,j}.incidentedAngle=atand(swath{i,j}.across_./abs(swath{i,j}.depth));
        end
    end
%------------   Simultaneous correction of the two main errors
    lastSoundSpeed=1501;%sound speed of last layer
    rollCorret=[1 -1 1 -1];%roll bias
%------------Processing of overlapping area in adjacent swaths
    figure("Position",[1957,150,1000,400],"Color",[1 1 1]);
    
    for i=1:1:col
        disp(i);
        clf('reset');
        plot(Xcoor,raster_real(:,i),'g.-');
        terrainProfile=raster_error(:,i);
        terrainProfile_corret  = terrainProfile;
        title([num2str(i) 'th']);
        hold on;
        
        if minDepthRaster<0
            minDepthRaster_positive = -minDepthRaster;
        end
        plot(Xcoor,raster_error(:,i),'r.-');
        plot( Xcoor(SwathFlag(i,:)),terrainProfile(SwathFlag(i,:)),'ms','MarkerSize',15);

        for j =1:swathCount
            swath{i,j}.corretIncidentedAngle=swath{i,j}.incidentedAngle+rollCorret(j);
            beamCount=length(swath{i,j}.depth);
            for k=1:1:beamCount
                theta=abs(swath{i,j}.corretIncidentedAngle(k));
                p=sind(theta) / surfaceSoundSpeed;
                acrossK= (sqrt(abs(1.0 - p*surfaceSoundSpeed * p*surfaceSoundSpeed)) - sqrt(abs(1.0 - p*lastSoundSpeed * p*lastSoundSpeed))) / (p*(lastSoundSpeed - surfaceSoundSpeed) / minDepthRaster_positive);
                timeK=(asin(p*lastSoundSpeed) - asin(p*surfaceSoundSpeed))*log(lastSoundSpeed / surfaceSoundSpeed)...
                    /...
                    (p*(lastSoundSpeed - surfaceSoundSpeed)*(lastSoundSpeed - surfaceSoundSpeed) / minDepthRaster_positive);
                acrossK=acrossK+(swath{i,j}.travelTime(k)-timeK)*lastSoundSpeed*sind(theta);%roll and refraction error correction
                depthK=minDepthRaster_positive+(swath{i,j}.travelTime(k)-timeK)*lastSoundSpeed*cosd(theta);
                if swath{i,j}.depth(k) < 0
                    swath{i,j}.corretDepth_(k)=-depthK;
                end
                if swath{i,j}.corretIncidentedAngle(k)<0
                    swath{i,j}.corretAcross_(k)=-acrossK;
                else
                    swath{i,j}.corretAcross_(k)=acrossK;
                end
            end%for k=1:1:beamCount
            % plot(swath{i,j}.x,swath{i,j}.corretDepth_,'b.-');
            swath{i,j}.corretDepth_Interp=interp1(swath{i,j}.corretAcross_,swath{i,j}.corretDepth_,swath{i,j}.across_,"linear","extrap");  
            
            if j>1
                terrainProfile_corret(swath{i,j}.startIndex+1:swath{i,j}.endIndex) = swath{i,j}.corretDepth_Interp(2:end);
            else
                terrainProfile_corret(swath{i,j}.startIndex:swath{i,j}.endIndex) = swath{i,j}.corretDepth_Interp;
            end
        end
        for j =2:swathCount
            if j==2
                flag1=2;
                flag2=3;
                flag3=4;
                flag4=5;
            elseif j==3
                flag1=6;
                flag2=7;
                flag3=8;
                flag4=9;
            elseif j==4
                flag1=10;
                flag2=11;
                flag3=12;
                flag4=13;
            end
            overLap_ascentFlag=SwathFlag(i,flag1);
            overLap_balanceStartFlag=SwathFlag(i,flag2);
            if overLap_balanceStartFlag <= overLap_ascentFlag
                overLap_balanceStartFlag = overLap_ascentFlag+3;
            end

            overLap_balanceEndFlag=SwathFlag(i,flag3);
            overLap_descentFlag=SwathFlag(i,flag4);
            if overLap_descentFlag <= overLap_balanceEndFlag
                overLap_descentFlag = overLap_balanceEndFlag+3;
            end
            overLap_terrain = terrainProfile(overLap_balanceStartFlag:overLap_balanceEndFlag);
            %______________balance stage
            LeftSwathDepth=swath{i,j-1}.corretDepth_Interp(overLap_ascentFlag-swath{i,j-1}.startIndex+1);
            RightSwathDepth=swath{i,j}.corretDepth_Interp(overLap_descentFlag-swath{i,j}.startIndex+1);
            depth_interp = linspace(LeftSwathDepth, RightSwathDepth, overLap_descentFlag-overLap_ascentFlag + 1)';
            depth_interp_balance = depth_interp(overLap_balanceStartFlag-overLap_ascentFlag+1 ...
                :end-(overLap_descentFlag-overLap_balanceEndFlag));
            % plot(Xcoor(overLap_balanceStartFlag:overLap_balanceEndFlag),depth_interp_balance,'k.-');
            maxDis = min(depth_interp_balance-overLap_terrain);
            minDis = max(depth_interp_balance-overLap_terrain);
            overLap_terrainCorrect=overLap_terrain+0.5*(maxDis+minDis);
            % plot(Xcoor(overLap_balanceStartFlag:overLap_balanceEndFlag),overLap_terrainCorrect,'.-');
            terrainProfile_corret(overLap_balanceStartFlag:overLap_balanceEndFlag) = overLap_terrainCorrect;
            %______________ascent stage
            if overLap_balanceStartFlag - overLap_ascentFlag > 1
                LeftSwathAscentIndex = overLap_ascentFlag:overLap_balanceStartFlag;%上升阶段的序列
                % LeftSwathAscentDepth = swath{i,j-1}.corretDepth_Interp(LeftSwathAscentIndex-swath{i,j-1}.startIndex+1);
                LeftSwathAscentDepth = terrainProfile (LeftSwathAscentIndex);
                startCorrect = terrainProfile_corret(overLap_ascentFlag);
                endDepth_corret = overLap_terrainCorrect(1);
                LeftSwathAscentDepthCorrect = startCorrect+(LeftSwathAscentDepth - LeftSwathAscentDepth(1))...
                                             /(LeftSwathAscentDepth(end)-LeftSwathAscentDepth(1))...
                                             *(endDepth_corret-startCorrect);
                terrainProfile_corret(overLap_ascentFlag+1:overLap_balanceStartFlag-1) = LeftSwathAscentDepthCorrect(2:end-1);
                % plot(Xcoor,terrainProfile_corret,'m.-');
            end
            %______________descent stage
            if overLap_descentFlag - overLap_balanceEndFlag > 1
                RightSwathAscentIndex = overLap_balanceEndFlag:overLap_descentFlag;
                RightSwathAscentDepth = terrainProfile (RightSwathAscentIndex);
                startCorrect = overLap_terrainCorrect(end);
                endDepth_corret = terrainProfile_corret(overLap_descentFlag);
                RightSwathAscentDepthCorrect = startCorrect+(RightSwathAscentDepth - RightSwathAscentDepth(1))...
                             /(RightSwathAscentDepth(end)-RightSwathAscentDepth(1))...
                             *(endDepth_corret-startCorrect);
                terrainProfile_corret(overLap_balanceEndFlag+1:overLap_descentFlag-1) = RightSwathAscentDepthCorrect(2:end-1);
                % plot(Xcoor,terrainProfile_corret,'k.-');
            end


        end
        plot(Xcoor,terrainProfile_corret,'k.-');
        
        box on;
        grid on;
        legend({'Befor errors correction','swath split point','Afer errors correction ','Real seafloor'});
        hold off;
        raster_correct(:,i) = terrainProfile_corret;
    end
save 'CorrectedSeafoor.mat' raster_correct raster_R;
end%function
