clc; clear all; close all;

% Change the current folder to the folder of this m-file.
if(~isdeployed)
    cd(fileparts(which(mfilename)));
end

%% 1st load data
idx = 8;

datadir = 'C:/Users/xiahaa/Documents/DTU/projectbank/image-stitching/data/fusion';

imu_prefix = strcat(datadir,'/imuRaw');
vicon_prefix = strcat(datadir,'/viconRot');
cam_prefix = strcat(datadir,'/cam');

    [imu_ts, imu_vals, vicon_ts, vicon_euler] = load_data(idx, imu_prefix, vicon_prefix);
    
if 0
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
imu_idx = 0;

Rvision = zeros(3,3,size(imgs,4));

imu_idx = findnn(imu_ts, imu_idx, img_ts(1));
Rvision(:,:,1) = euler_to_rot(gyro_euler(3,imu_idx),gyro_euler(2,imu_idx),gyro_euler(1,imu_idx));

img1 = imgs(:,:,:,1);
img1d = single(rgb2gray(img1));

Homos = zeros(3,3,size(imgs,4)-1);
flags = zeros(1,size(imgs,4));

doHomo = 0;

if doHomo == 0 || doHomo == 1
    load('homo8.mat','Homos','flags');
end

wfov = pi / 3;
hfov = pi / 4;
width = size(imgs,2);
height = size(imgs,1);
fx = width * 0.5 / tan(wfov/2);
fy = height * 0.5 / tan(hfov/2);

K = [fx 0 width*0.5;0 fy height*0.5;0 0 1];

Ric = [0 0 1;-1 0 0;0 -1 0];

figure
wrong = 0;
for i = 2:size(imgs,4)
    img2 = imgs(:,:,:,i);    
    img2d = single(rgb2gray(img2));

    %% last imu
    prevImu = euler_to_rot(gyro_euler(3,imu_idx),gyro_euler(2,imu_idx),gyro_euler(1,imu_idx));
    prevTime = imu_ts(imu_idx);
    %% new imu_index
    imu_idx = findnn(imu_ts, imu_idx, img_ts(i));
    if imu_idx >= length(imu_ts)
        break;
    end
    
    %% new imu
    Rimu = euler_to_rot(gyro_euler(3,imu_idx),gyro_euler(2,imu_idx),gyro_euler(1,imu_idx));
    imuTime = imu_ts(imu_idx);
    
    %% if do homography estimation
    if doHomo == 1 || doHomo == 2
        [feature1,descriptor1] = vl_sift(img1d) ;
        [feature2,descriptor2] = vl_sift(img2d) ;

        ds2 = single(descriptor2);
        ds1 = single(descriptor1);
        kdtree = vl_kdtreebuild(ds2);
        matches = zeros(size(descriptor1,2),2);
        k = 1;
        thresh = 2;
        for ii = 1:size(descriptor1,2)
            [index, distance] = vl_kdtreequery(kdtree, ds2, ds1(:,ii), ...
                                               'NumNeighbors', 2, 'MaxComparisons', 15);
            if distance(1)*thresh > distance(2) continue; end
            matches(k,:) = [ii,index(1)];
            k = k + 1;
        end
        matches(k:end,:) = [];
        matches = matches';
        
        x1h = tohomogeneous(feature1(1:2,matches(1,:)));
        x2h = tohomogeneous(feature2(1:2,matches(2,:)));

        if doHomo == 2
            ransac.pinlier = 0.99;
            ransac.estt_fun = @HestWithNormalization;%plane_estimation
            ransac.eval_fun = @reprojectionError;%dist2plane
            ransac.maxiter = 1e3;
            ransac.threshold = 6;
            ransac.inliers = [];
            ransac.minimumset = 4;

            if size(x1h,2) < ransac.minimumset
                img1 = img2;
                img1d = single(rgb2gray(img1));
                continue;
            end

            result = ransac_routine_homo(x1h, x2h, ransac);
            Homos(:,:,i) = result.params;
            flags(i) = 1;
            x1h = x1h(:,result.inliers);
            x2h = x2h(:,result.inliers);
            
            imshow1 = cat(2, img1, img2);
            imshow(imshow1);hold on;

            plot(x1h(1,:),x1h(2,:), 'ro','MarkerSize',5);
            plot(x2h(1,:)+size(img1,2),x2h(2,:), 'bo','MarkerSize',5);

            shift = size(img1,2);
            cmap = jet(32);
            k = 1;
            for j = 1:size(x1h,2)
                ptdraw = [x1h(2,j), x1h(1,j);
                          x2h(2,j), x2h(1,j)+shift];
                plot(ptdraw(:,2),ptdraw(:,1),'LineStyle','-','LineWidth',0.5,'Color',cmap(k,:));
                k = mod(k+1,32);if k == 0 k = 1;end
            end
            pause(0.1);
        else
            imshow1 = cat(2, img1, img2);
            imshow(imshow1);hold on;

            plot(x1h(1,:),x1h(2,:), 'ro','MarkerSize',5);
            plot(x2h(1,:)+size(img1,2),x2h(2,:), 'bo','MarkerSize',5);
            if flags(i) == 1
                x3h = Homos(:,:,i) * x1h;
                x3h = x3h./x3h(3,:);
                plot(x3h(1,:)+size(img1,2),x3h(2,:), 'g+','MarkerSize',5);
            end
            pause(0.1);
        end
    else
        %% only do homography decomposition 
        if flags(i) == 1
            H1 = Homos(:,:,i);
            H1 = H1./H1(3,3);
            H1 = normalization_homography(H1);
            R2 = [];
            [R2,t2,n2] = homo_decom_malis(H1);
            
            % delta rotation from last camera to current camera measured by
            % imu
            dRimu = Ric'*Rimu'*prevImu*Ric;
            
            % find the most correct one
            minval = 1e6;
            minval = rot_err(dRimu,R2(:,:,1));
            minid = 1;
            minR = R2(:,:,1);
            for ii = 2:size(R2,3)
                err = rot_err(dRimu,R2(:,:,ii));
                if err < minval
                    minval = err;
                    minid = ii;
                    minR = R2(:,:,ii);
                end
            end
%             disp(minval);
%             disp(minid)
%             minvals(i) = minval;
            if minval > 0.1
                wrong = wrong + 1;
                disp(wrong);
                Rvision(:,:,i) = Rimu;
                minRvs(:,:,i) = dRimu;
                minRv_ts(i) = img_ts(i);
            else
%             Rvision(:,:,i) = prevImu;
%                 minRvs(:,:,i) = minR;
%                 minRv_ts(i) = img_ts(i);
%                 minRis(:,:,i) = dRimu;
                dR = fusionVisionInertial(minR, dRimu, img_ts(i)-img_ts(i-1));

                Rvision(:,:,i) = Rvision(:,:,i-1) * (Ric * dR' * Ric');%Rvision(:,:,i-1) * (Ric * R2(:,:,minid)' * Ric')';
            end
        else
            Rvision(:,:,i) = Rimu;
        end
    end
    img1 = img2;
    img1d = single(rgb2gray(img1));
end

save('vision_meas.mat','Rvision');

% n = size(minRvs,3);
% homo_euler = zeros(3,n);
% homo_euler2 = zeros(3,n);
% for i = 1:1:n
%     [r,p,y] = rot_to_euler(minRvs(:,:,i));
%     homo_euler(:,i) = [r;p;y];
%     [r,p,y] = rot_to_euler(minRis(:,:,i));
%     homo_euler2(:,i) = [r;p;y];
% end
% figure;
% subplot(3,1,1);
% plot(homo_euler(1,:));hold on;grid on;
% plot(homo_euler2(1,:));
% legend({'SO3','Vision'},'FontName','Arial','FontSize',10);
% title('Attitude Estimation','FontName','Arial','FontSize',20);grid on;
% subplot(3,1,2);plot(homo_euler(2,:));hold on;
% plot(homo_euler2(2,:));grid on;
% subplot(3,1,3);plot(homo_euler(3,:));hold on;
% plot(homo_euler2(3,:));grid on;



n = size(Rvision,3);
vision_euler = zeros(3,n);
for i = 1:1:n
    [r,p,y] = rot_to_euler(Rvision(:,:,i));
    vision_euler(:,i) = [r;p;y];
end
figure;
subplot(3,1,1);plot(imu_ts, gyro_euler(1,:),'r-');hold on;
plot(img_ts, vision_euler(1,:),'g-');
plot(vicon_ts, vicon_euler(1,:),'b-');
legend({'SO3','Vision','VICON'},'FontName','Arial','FontSize',10);
title('Attitude Estimation','FontName','Arial','FontSize',20);grid on;
subplot(3,1,2);plot(imu_ts, gyro_euler(2,:),'r-');hold on;plot(img_ts, vision_euler(2,:),'g-');grid on;
plot(vicon_ts, vicon_euler(2,:),'b-');
subplot(3,1,3);plot(imu_ts, gyro_euler(3,:),'r-');hold on;plot(img_ts, vision_euler(3,:),'g-');grid on;
plot(vicon_ts, vicon_euler(3,:),'b-');


if doHomo == 1
    save('homo.mat','Homos','flags');
end


function dR = fusionVisionInertial(dRv, dRi, dT)
    so3i = logm(dRv'*dRi);
    wi = ([-so3i(2,3);so3i(1,3);-so3i(1,2)]);
    wi = wi.*[-1;1;-1];
    fskew = @(x) ([0 -x(3) x(2);x(3) 0 -x(1);-x(2) x(1) 0]);
    kp = 1.5;
    dt = dT;
    so3 = fskew((wi * kp));
    dR = dRi*expm(so3.*dt);%% dt: integration time
    RtR = (dR)'*(dR);
    E = RtR - eye(3);
    if max(abs(E)) > 1e-6
        %% orthogonization
        [U, ~, V] = svd(dR); 
        dR = U*V';
        if det(dR)<0, dR = U*diag([1 1 -1])*V'; end 
    end
end

function err = rot_err(R1,R2)
    e1 = logm(R2*R1');
    err = norm([-e1(2,3);e1(1,3);-e1(1,2)]);
end

function imu_idx = findnn(imu_ts, imu_idx, img_ts)
    if imu_idx >= length(imu_ts)
        return;
    else
        timediff = abs(img_ts - imu_ts(imu_idx+1:end));
        [~,minid] = min(timediff);
    end
    % get nearest rotation matrix
    imu_idx = imu_idx + minid;
end

function Hn = normalization_homography(H)
    [~,S,~]=svd(H);
    Hn = H./S(2,2);
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

function varargout = homo_decom_malis(varargin)
%% decompose homography matrix using olivier algo.
    H = varargin{1};
    format long;
    
    S = H'*H-eye(3);
    
    s11 = S(1,1);s12 = S(1,2);s13 = S(1,3);s22 = S(2,2);s23 = S(2,3);s33 = S(3,3);
%     s21 = S(2,1);
%     s31 = S(3,1);s32 = S(3,2);
    
    Ms11 = s23*s23 - s22*s33;
    Ms22 = s13*s13 - s11*s33;
    Ms33 = s12*s12 - s11*s22;
    
    Ms12 = s23*s13 - s12*s33;
    Ms13 = s22*s13 - s12*s23;
    Ms23 = s12*s13 - s11*s23;
    
    vscalar = 2*sqrt(1+trace(S)-Ms11-Ms22-Ms33);
    tenorm = sqrt(2+ trace(S)-vscalar);

    %% intermediate varibales
    if abs(s22) > 1e-6
        epsilon13 = signc(Ms13);
        nea = [s12+sqrt(Ms33);s22;s23-epsilon13*sqrt(Ms11)];
        neb = [s12-sqrt(Ms33);s22;s23+epsilon13*sqrt(Ms11)];
        
        na = nea./norm(nea);
        nb = neb./norm(neb);
        
        epsilons22 = signc(s22);
        rho = sqrt(tenorm*tenorm+2*vscalar);
        tas = tenorm*0.5*(epsilons22*rho*nb-tenorm*na);
        tbs = tenorm*0.5*(epsilons22*rho*na-tenorm*nb);
        
    elseif abs(s11) > 1e-6
        epsilon23 = signc(Ms23);
        nea = [s11;s12+sqrt(Ms33);s13+epsilon23*sqrt(Ms22)];
        neb = [s11;s12-sqrt(Ms33);s13-epsilon23*sqrt(Ms22)];
        
        na = nea./norm(nea);
        nb = neb./norm(neb);
        
        epsilons11 = signc(s11);
        rho = sqrt(tenorm*tenorm+2*vscalar);
        tas = tenorm*0.5*(epsilons11*rho*nb-tenorm*na);
        tbs = tenorm*0.5*(epsilons11*rho*na-tenorm*nb);
        
    elseif abs(s33) > 1e-6
        epsilon12 = signc(Ms12);
        nea = [s13+epsilon12*sqrt(Ms22);s23+sqrt(Ms11);s33];
        neb = [s13-epsilon12*sqrt(Ms22);s23-sqrt(Ms11);s33];
        
        na = nea./norm(nea);
        nb = neb./norm(neb);
        
        epsilons33 = signc(s33);
        rho = sqrt(tenorm*tenorm+2*vscalar);
        tas = tenorm*0.5*(epsilons33*rho*nb-tenorm*na);
        tbs = tenorm*0.5*(epsilons33*rho*na-tenorm*nb);
    else
        error('no solution');
    end
    
    
    n = [na nb];
    R(:,:,1) = H*(eye(3)-2/vscalar*tas*na');
    R(:,:,2) = H*(eye(3)-2/vscalar*tbs*nb');
    
    R(:,:,1) = projectOnSO3(R(:,:,1));
    R(:,:,2) = projectOnSO3(R(:,:,2));
    
    t(:,1) = R(:,:,1)*tas;
    t(:,2) = R(:,:,2)*tbs;
    
    n = [n -na -nb];
    R(:,:,3) = R(:,:,1);
    R(:,:,4) = R(:,:,2);
    
    t(:,3) = -t(:,1);
    t(:,4) = -t(:,2);
    
    valid1 = ones(1,size(t,2));

    if nargin == 2
        m1 = varargin{2};
        for i = 1:size(t,2)
            if n(:,i)'*m1(:,1) > 0
                valid1(i) = 1;
            else
                valid1(i) = 0;
            end
        end
    end
    valid1 = valid1 == 1;
    Rf = R(:,:,valid1);
    tf = t(:,valid1);
    nf = n(:,valid1);
%         za1 = s12+sqrt(Ms33)/s22; zb1 = s12-sqrt(Ms33)/s22;
%         za3 = s23-epsilon13*sqrt(Ms11)/s22;zb3 = s23+epsilon13*sqrt(Ms11)/s22;
%         aa = 1 + za1*za1+za3*za3;
%         ab = 1 + zb1*zb1+zb3*zb3;
%         b = 2+trace(S);    
    
    varargout{1} = Rf;
    varargout{2} = tf;
    varargout{3} = nf;
end

function R = projectOnSO3(R)
    RtR = (R)'*(R);
    E = RtR - eye(3);
    if max(abs(E)) > 1e-6
        %% orthogonization
        [U, ~, V] = svd(R); 
        R = U*V';
        if det(R)<0, R = U*diag([1 1 -1])*V'; end 
    end
end

function s = signc(a)
    if a > 0 || abs(a) < 1e-6
        s = 1;
    else
        s = -1;
    end
end