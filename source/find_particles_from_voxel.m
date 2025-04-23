function TotPosArr = find_particles_from_voxel(vox)
%% Main Polynomial Tracing
% code to perform polynomial tracing on reconstruction volume
% reconstruction volume should be upsampled by 3*3*3 linear interpolation


% each local maximum is fitted with 9*9*9 voxels (3*3*3 before interpolation)
% by 4th order polynomial equation to get the position

% % Th: intensity threshold for the local maxima pixel
% local maxima with intensity less than this value will not be traced
% because they are way too weak to become actual atoms
MaxIter = 14;       
CritIter = 7;         
minDist  = 2;     
SearchRad = 3;
% initialize the parameters
Q = 0.5;     Alpha = 1;
cropHalfSize = SearchRad;


% pre-search particles;
pos = getfallingedge(vox);
fprintf('Falling edge search done.\n')
maxXYZ = vox2pos(pos);

t = maxXYZ(:,1) <= cropHalfSize | maxXYZ(:,1) >= size(vox,1)-cropHalfSize;
t = t | maxXYZ(:,2) <= cropHalfSize | maxXYZ(:,2) >= size(vox,2)-cropHalfSize;
t = t | maxXYZ(:,3) <= cropHalfSize | maxXYZ(:,3) >= size(vox,3)-cropHalfSize;
maxXYZ(t, :) = [];

% get polynomial power array
fitCoeff = [];
for i=0:4
    for j=0:4
        for k=0:4
            if i+j+k <= 4                
                if max([i j k]) == 4
                    fitCoeff(end+1,:) = [i j k -1];                %#ok<SAGROW>
                else                
                    fitCoeff(end+1,:) = [i j k 0];                 %#ok<SAGROW>
                end
            end
        end
    end
end


[X,Y,Z]      = ndgrid(-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize);
SphereInd    = find(X.^2+Y.^2+Z.^2 <=(SearchRad+0.5)^2);
XYZdata.X    = X(SphereInd); XYZdata.Y = Y(SphereInd); XYZdata.Z = Z(SphereInd);

Orders       = fitCoeff(:,1:3);
PosArr       = zeros(size(maxXYZ));
TotPosArr    = zeros(size(maxXYZ));

exitFlagArr  = zeros(1, size(maxXYZ,1));
CoeffArr     = repmat(fitCoeff(:,4),[1 size(maxXYZ,1)]);


% perform the main tracing loop
for i=1:size(maxXYZ,1)
    endFlag = 0;
    consecAccum = 0;
    iterNum = 0;
    while ~endFlag    
        iterNum = iterNum + 1;
        if iterNum>MaxIter
          exitFlagArr(i) = -4;
          endFlag = 1;
        end
        cropXind = maxXYZ(i,1) + (-cropHalfSize:cropHalfSize);
        cropYind = maxXYZ(i,2) + (-cropHalfSize:cropHalfSize);
        cropZind = maxXYZ(i,3) + (-cropHalfSize:cropHalfSize);

        cropVol = vox(cropXind,cropYind,cropZind);

        Pos = PosArr(i,:);
        ydata = cropVol(SphereInd);
        GaussWeight = exp(-1*Alpha*( (X(SphereInd)-Pos(1)).^2 + (Y(SphereInd)-Pos(2)).^2 + (Z(SphereInd)-Pos(3)).^2 ) / cropHalfSize^2 );
        
        %fun = @(p,xdata) calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z,Pos,Orders,p).*GaussWeight;

        %opts = optimset('Display','off');
        
        %[p1,fminres1] = lsqcurvefit(fun,CoeffArr(:,i),XYZdata,cropVol(SphereInd).*GaussWeight,[],[],opts);
        options = optimset('fminsearch');
        options = optimset(options, 'TolX',0.1E-8);
        %        options = optimset(options, 'PlotFcns',@optimplotx);
        options = optimset(options, 'MaxIter',100);
        options = optimset(options, 'MaxFunEvals', 100);
        
        p1 = fminsearch(@(x) myfitfunc(x, ydata.*GaussWeight, Pos, XYZdata, Orders, GaussWeight),CoeffArr(:,i), options);

        CoeffArr(:,i) = p1;
        
        [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,CoeffArr(:,i));
        if dX ==-100 && dY == -100 && dZ == -100
            exitFlagArr(i) = -1;
            endFlag = 1;
        else
            maxedShift = max([dX dY dZ],-1*[Q Q Q]);
            minedShift = min(maxedShift,[Q Q Q]);
            PosArr(i,:) = PosArr(i,:) + minedShift;
            if max(abs(PosArr(i,:))) > cropHalfSize
                exitFlagArr(i) = -2;
                endFlag = 1;
            elseif max(abs(minedShift)) < Q
                if consecAccum == CritIter-1
                    goodAtomTotPos = TotPosArr(1:i-1,:);
                    goodAtomTotPos = goodAtomTotPos(exitFlagArr(1:i-1)==0,:);
                    Dist = sqrt(sum((goodAtomTotPos - repmat(PosArr(i,:)+maxXYZ(i,:),[size(goodAtomTotPos,1) 1])).^2,2));
                    if min(Dist) < minDist
                      exitFlagArr(i) = -3;
                    else
                      TotPosArr(i,:) = PosArr(i,:) + maxXYZ(i,:);
                    end
                    endFlag = 1;
                else
                    consecAccum = consecAccum + 1;
                end
            else
                consecAccum = 0;
            end
        end
    end
    fprintf(1,'peak %d, flag %d \n',i,exitFlagArr(i));
end

function err = myfitfunc(coeff, data, p, xdata, Orders, GaussWeight)
y = calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z, p, Orders, coeff).*GaussWeight;
err = sum((y-data).^2)/numel(p);
