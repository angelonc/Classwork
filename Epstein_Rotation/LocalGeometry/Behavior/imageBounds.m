function imageBounds

imDir = ['~/Documents/Classwork/Epstein_Rotation/LocalGeometry/' ...
         'Maze_Unity/Screenshots/'];
fn = '~/Documents/Classwork/Epstein_Rotation/LocalGeometry/Behavior/Configs/cropped_ims.mat';
im_files = dir([imDir '*.png']);

plot = 1;
crop = 1;

% If file not made
if ~exist(fn,'file')
    % For each image
    rect = zeros(length(im_files),4);
    for i = 1:length(im_files)
        clear im tmp
        
        % Set to grayscale
        tmp = imread([imDir im_files(i).name]);
        im = mean(tmp,3);
        % imagesc(im);
        
        % Find bounding box
        [r,c] = find(diff(im) > 0);
        
        rect(i,1) = min(c);
        rect(i,2) = min(r);
        rect(i,3) = max(c);
        rect(i,4) = max(r);
        
        if plot
            im(:,rect(i,1):rect(i,1)+2) = 0;
            im(rect(i,2):rect(i,2)+2,:) = 0;
            im(:,rect(i,3):rect(i,3)+2) = 0;
            im(rect(i,4):rect(i,4)+2,:) = 0;
            imagesc(im);
            pause(.0001);
        end
        
        % Find value of grey background
        if i == 1
            bg = mode(im(:));
        end
        
        % Crop
        if crop
            cRect = [rect(i,1) rect(i,2) rect(i,3)-rect(i,1) ...
                     rect(i,4)-rect(i,2)];
            ims{i} = imcrop(tmp,cRect);
        end
        fprintf('%03d\n',i);
    end
    
    save(fn,'im_files','rect','bg','ims');
    disp('Done!');
end

    






    