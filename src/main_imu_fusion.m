clc; clear all; close all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
    cd(fileparts(which(mfilename)));
end

%% 1st load data
idx = 10;

datadir = 'C:/Users/xiahaa/Documents/DTU/projectbank/image-stitching/data/fusion';

% imu_prefix = strcat(datadir,'/imu/imuRaw');
% vicon_prefix = '/vicon/viconRot';
% cam_prefix = "cam/cam"
imu_prefix = strcat(datadir,'/imuRaw');
vicon_prefix = strcat(datadir,'/viconRot');


[imu_ts, imu_vals, vicon_vs, vicon_euler] = load_data(idx, imu_prefix, vicon_prefix);

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
    acc_estimate(:,i) = R' * a_ref;
end

figure;
subplot(3,1,1);plot(imu_ts, gyro_euler(1,:));hold on;plot(vicon_vs, vicon_euler(1,:));
subplot(3,1,2);plot(imu_ts, gyro_euler(2,:));hold on;plot(vicon_vs, vicon_euler(2,:));
subplot(3,1,3);plot(imu_ts, gyro_euler(3,:));hold on;plot(vicon_vs, vicon_euler(3,:));
legend({'SO3','VICON'})

figure;
subplot(3,1,1);plot(imu_ts, acc_truth(1,:));hold on;plot(imu_ts, acc_estimate(1,:));
subplot(3,1,2);plot(imu_ts, acc_truth(2,:));hold on;plot(imu_ts, acc_estimate(2,:));
subplot(3,1,3);plot(imu_ts, acc_truth(3,:));hold on;plot(imu_ts, acc_estimate(3,:));


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
    vicon =  load(strcat(vicon_prefix,num2str(idx),'.mat'));
    vicon_vals = vicon.rots;
    vicon_ts = vicon.ts;
    n = size(vicon_vals,3);
    vicon_euler = zeros(3,n);
    for i = 1:1:n
        [r,p,y] = rot_to_euler(vicon_vals(:,:,i));
        vicon_euler(:,i) = [r;p;y];
    end
    
end

function [y,p,r] = rot_to_euler(rot)

    p =  sin(rot(1,3));
    r = -atan2(rot(2,3),rot(3,3));
    y = -atan2(rot(1,2),rot(1,1));

end