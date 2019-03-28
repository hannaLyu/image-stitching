clc; clear all; close all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
    cd(fileparts(which(mfilename)));
end

%% 1st load data
idx = 11;

datadir = 'C:/Users/xiahaa/Documents/DTU/projectbank/image-stitching/data/fusion';

% imu_prefix = strcat(datadir,'/imu/imuRaw');
% vicon_prefix = '/vicon/viconRot';
% cam_prefix = "cam/cam"
imu_prefix = strcat(datadir,'/imuRaw');
vicon_prefix = strcat(datadir,'/viconRot');
cam_prefix = strcat(datadir,'/cam');

    [imu_ts, imu_vals, vicon_ts, vicon_euler] = load_data(idx, imu_prefix, vicon_prefix);
if 1
    n = size(imu_ts,2);
    gyro_euler = zeros(3,n);  % represent orientation in euler angles
    acc_estimate = zeros(3,n); % represent orientation in euler angles
    acc_truth = imu_vals(1:3,:);

    R = eye(3);
    kp = 10;
    wb = 0;

    a_ref = [0;0;1];
    for  i = 1:n
        % extract sensor data
        acc = imu_vals(1:3,i);
        gyro = imu_vals(4:6,i);

        if i == n
            dt = mean(imu_ts(end-10:end) - imu_ts(end-11:end-1));
        else
            dt = imu_ts(i+1) - imu_ts(i);
        end
        %% nomialal
        an = R' * a_ref;
        wmeas = cross(acc, an); %% am mm: measurement
        %% Explicit complementary filter on SO3
        fskew = @(x) ([0 -x(3) x(2);x(3) 0 -x(1);-x(2) x(1) 0]);
        so3 = fskew((gyro - wb + wmeas * kp));
        R = R*expm(so3.*dt);%% dt: integration time
        RtR = (R)'*(R);
        E = RtR - eye(3);
        if max(abs(E)) > 1e-6
            %% orthogonization
            [U, ~, V] = svd(R); 
            R = U*V';
            if det(R)<0, R = U*diag([1 1 -1])*V'; end 
        end
        [r,p,y] = rot_to_euler(R);
        gyro_euler(:,i) = [r;p;y];
        R1 = euler_to_rot(gyro_euler(3,i),gyro_euler(2,i),gyro_euler(1,i));
        e1 = logm(R*R1');
        err(i) = norm([-e1(2,3);e1(1,3);-e1(1,2)]);
        acc_estimate(:,i) = R' * a_ref;
    end
    save('process1.mat','gyro_euler','acc_estimate','acc_truth');
else
    load('process1.mat');
end
if ~isempty(vicon_ts)
    figure;
    subplot(3,1,1);plot(imu_ts-imu_ts(1), gyro_euler(1,:));hold on;
    plot(vicon_ts-vicon_ts(1), vicon_euler(1,:));legend({'SO3','VICON'},'FontName','Arial','FontSize',10);
    title('Attitude Estimation','FontName','Arial','FontSize',20);grid on;
    subplot(3,1,2);plot(imu_ts-imu_ts(1), gyro_euler(2,:));hold on;plot(vicon_ts-vicon_ts(1), vicon_euler(2,:));grid on;
    subplot(3,1,3);plot(imu_ts-imu_ts(1), gyro_euler(3,:));hold on;plot(vicon_ts-vicon_ts(1), vicon_euler(3,:));grid on;
end
figure;
subplot(3,1,1);plot(imu_ts-imu_ts(1), acc_truth(1,:));hold on;plot(imu_ts-imu_ts(1), acc_estimate(1,:));grid on;
legend({'Truth','Prediction'},'FontName','Arial','FontSize',10);title('Measurement Comparison','FontName','Arial','FontSize',20);
subplot(3,1,2);plot(imu_ts-imu_ts(1), acc_truth(2,:));hold on;plot(imu_ts-imu_ts(1), acc_estimate(2,:));grid on;
subplot(3,1,3);plot(imu_ts-imu_ts(1), acc_truth(3,:));hold on;plot(imu_ts-imu_ts(1), acc_estimate(3,:));grid on;

%% load cam
cam = load(strcat(cam_prefix,num2str(idx),'.mat'));
imgs = cam.cam;
img_ts = cam.ts;

%% rough estimation of the intrinsics
wfov = pi / 3;
hfov = pi / 4;
width = size(imgs,2);
height = size(imgs,1);
fx = width * 0.5 / tan(wfov/2);
fy = height * 0.5 / tan(hfov/2);
f = (fx+fy)*0.5;

imu_idx = 0;

[yy,xx] = meshgrid(1:height,1:width);
xx = vec(xx');
yy = vec(yy');
xx = (xx - round(width*0.5));
yy = (yy - round(height*0.5));
zz = f.*ones(numel(xx),1);
[cx,cy,cz] = cartesian_to_cylindrical(xx,yy,zz);

ncy = round(pi / wfov * f);
ncx = round(pi*0.5 / hfov * f);

frame = zeros(round(pi / hfov * f)+1, round(pi*2 / wfov * f)+1, size(imgs,3));
weights = zeros(size(frame,1),size(frame,2));

weight = zeros(height,width);
weight(1,:) = 1;weight(end,:) = 1;weight(:,1) = 1;weight(:,end) = 1;
weight = bwdist(weight,'euclidean');%maxdist = max(dist1(:));dist1 = (maxdist+1) - dist1;
Ric = [0 -1 0;0 0 -1;1 0 0];
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
for i = 1:size(imgs,4)
    img = imgs(:,:,:,i);
    imgt = img_ts(i);
    
    if imu_idx >= length(imu_ts)
        break;
    else
        timediff = abs(img_ts(i) - imu_ts(imu_idx+1:end));
        [minval,minid] = min(timediff);
%         disp(minval);
    end
    % get nearest rotation matrix
    imu_idx = imu_idx + minid;
    R1 = euler_to_rot(gyro_euler(3,imu_idx),gyro_euler(2,imu_idx),gyro_euler(1,imu_idx));
    
    
    % transformation
    rotated_cylind = R1 * [cx';cy';cz';];
    
    azimuth = -atan2(rotated_cylind(2,:),rotated_cylind(1,:));
    altitude = -rotated_cylind(3,:);%-atan2(rotated_cylind(3,:),sqrt(rotated_cylind(1,:).^2+rotated_cylind(2,:).^2));
    
    projx = round(azimuth ./ wfov .* f + ncy)' + 1;
    projy = round(altitude ./ hfov .* f + ncx)' + 1;
    
    for j = 1:size(img,3)
        [frame(:,:,j)] = merge(im2double(img(:,:,j)),frame(:,:,j),weight,projx,projy);
    end
%     weights((projx-1).*size(weights,1)+projy) = weights((projx-1).*size(weights,1)+projy) + weight(:);
%     
%     for j = 1:size(img,3)
%         display(:,:,j) = frame(:,:,j) ./ (weights + 1e-3);
%     end
%     weights = weights ./(weights + 1e-3);
    imshow(frame);
    
    % Write to the GIF File 
    % Capture the plot as an image 

%     imwrite(frame,strcat('../results/fusion/',num2str(i),'.png')); 

    
end

for j = 1:size(img,3)
    frame(:,:,j) = frame(:,:,j) ./ (weights + 1e-3);
end

    
function [dst]= merge(src,dst,weight,projx,projy)
    dst((projx-1).*size(dst,1)+projy) = src(:);%dst((projx-1).*size(dst,1)+projy) + weight(:).*
end

function [cx,cy,cz] = cartesian_to_cylindrical(x,y,z)
    r = sqrt(x.*x+z.*z);
    cx = z./r;
    cy = -x./r;
    cz = -y./r;
end

function [imu_ts, imu_vals, vicon_ts, vicon_euler] = load_data(idx, imu_prefix, vicon_prefix)
    imu = load(strcat(imu_prefix,num2str(idx),'.mat'));
    imu_vals = imu.vals;
    imu_ts = imu.ts;
    
    % scale and bias based on IMU reference
    acc_x = -imu_vals(1,:);
    acc_y = -imu_vals(2,:);
    acc_z =  imu_vals(3,:);
    acc = [acc_x;acc_y;acc_z];

    Vref = 3300;
    acc_sensitivity = 330;
    acc_scale_factor = Vref/1023/acc_sensitivity;
    acc_bias = mean(acc(:,1:100),2) - [0;0;1]./acc_scale_factor;
    acc = (acc-repmat(acc_bias,1,size(acc,2))).*acc_scale_factor;
    
    gyro_x = imu_vals(5,:);
    gyro_y = imu_vals(6,:);
    gyro_z = imu_vals(4,:);
    gyro = [gyro_x;gyro_y;gyro_z];
    gyro_sensitivity = 3.33;
    gyro_scale_factor = Vref/1023/gyro_sensitivity;
    gyro_bias = mean(gyro(:,1:100),2);
    gyro = (gyro-repmat(gyro_bias,1,size(gyro,2))).*gyro_scale_factor*(pi/180);

    imu_vals = [acc;gyro];
    
    % load vicon
    try
        vicon =  load(strcat(vicon_prefix,num2str(idx),'.mat'));
        vicon_vals = vicon.rots;
        vicon_ts = vicon.ts;
        n = size(vicon_vals,3);
        vicon_euler = zeros(3,n);
        for i = 1:1:n
            [r,p,y] = rot_to_euler(vicon_vals(:,:,i));
            vicon_euler(:,i) = [r;p;y];
        end
    catch
        vicon_ts = [];
        vicon_euler = [];
    end
end

function rot = euler_to_rot(roll, pitch, yaw)
    % from eular angle to rotation matrix, OK

    rot1 = [cos(yaw) sin(yaw) 0;-sin(yaw) cos(yaw) 0; 0 0 1];
    rot2 = [cos(pitch) 0 -sin(pitch);0 1 0; sin(pitch) 0 cos(pitch)];
    rot3 = [1 0 0;0 cos(roll) sin(roll);0 -sin(roll) cos(roll)];

    rot = rot3*rot2*rot1;
end


function [y,p,r] = rot_to_euler(rot)
    format long;
    p =  -asin(rot(1,3));
    r = atan2(rot(2,3),rot(3,3));
    y = atan2(rot(1,2),rot(1,1));

end