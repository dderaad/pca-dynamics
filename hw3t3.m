%% (test 3) Horizontal Displacement
clear variables; close all; clc;
%% Cam 1_3

load cam1_3.mat

datasize = size(vidFrames1_3); % height, width, color channels, length

x_a = zeros([datasize(4) 1]); y_a = x_a;

for j=1:datasize(4)
    gryfrm = rgb2gray(vidFrames1_3(:,:,:,j));
    frm = gryfrm>=245;
    frm(1:200,:) = 0; % remove artifact from top
    frm(:,1:250) = 0; % remove artifact from left side
    frm(:,end-100:end) = 0; % remove artifact from left side
    %imshow(frm, []); 
    %hold on
    [x_a(j), y_a(j)] = centroid(frm);
    %scatter(xcom, ycom, 50, 'x');
    %hold off
    %drawnow
end

clear vidFrames1_3;
'Cam1_3 done.'

%plot(y_a)

%% Cam 2_3

load cam2_3.mat;

x_b = zeros([datasize(4) 1]); y_b = x_b;

for j=1:datasize(4)
    gryfrm = rgb2gray(vidFrames2_3(:,:,:,j));
    frm = gryfrm>=250;
    frm(1:175,:) = 0; % remove artifact from top
    %imshow(frm, []); 
    %hold on
    [x_b(j), y_b(j)] = centroid(frm);
    %scatter(xcom, ycom, 50, 'x');
    %hold off
    drawnow
end

clear vidFrames2_1;
'Cam2_3 done.'

%plot(y_b)

%% Cam 3_3
load cam3_3.mat;

datasize = size(vidFrames3_3); % height, width, color channels, length

x_c = zeros([datasize(4) 1]); y_c = x_c;

for j=1:datasize(4)
    gryfrm = rgb2gray(vidFrames3_3(:,:,:,j));
    frm = gryfrm>=245;
    frm(:,1:200) = 0; % remove artifact from left side
    %imshow(frm, []); 
    %hold on
    [x_c(j), y_c(j)] = centroid(frm);
    %scatter(xcom, ycom, 50, 'x');
    %hold off
    drawnow
end

clear vidFrames3_1;
'Cam3_3 done.'

%plot(x_c)

%% (test 4) Horizontal Displacement Analysis

% Not all sensors started recording at the same time.
% We need to synchronize the recordings. This means we will have
% to synchronize the 1st shift in direction and trim frames not shared 
% by all videos.

t_a = 1:length(x_a);
t_b = 1:length(x_b);
t_c = 1:length(x_c);

figure(3)
subplot(3,1,1)
[pksa, loca] = findpeaks(y_a, 'MinPeakDistance', 35);
findpeaks(y_a, 'MinPeakDistance', 35)
subplot(3,1,2)
[pksb, locb] = findpeaks(y_b, 'MinPeakDistance', 35); % 2nd peak
findpeaks(y_b, 'MinPeakDistance', 35)
subplot(3,1,3)
[pksc, locc] = findpeaks(x_c, 'MinPeakDistance', 35);
findpeaks(x_c, 'MinPeakDistance', 35)

syncpeaks = [loca(2) locb(2) locc(2)];

chopinit = syncpeaks - min(syncpeaks); % 
achop = 1:chopinit(1);
x_a(achop) = []; y_a(achop) = []; t_a(achop) = [];
bchop = 1:chopinit(2);
x_b(bchop) = []; y_b(bchop) = []; t_b(bchop) = [];
cchop = 1:chopinit(3);
x_c(cchop) = []; y_c(cchop) = []; t_c(cchop) = [];

maxvidlength = min([length(t_a) length(t_b) length(t_c)]);
t = 1:maxvidlength;

x_a = x_a(t); x_b = x_b(t); x_c = x_c(t);
y_a = y_a(t); y_b = y_b(t); y_c = y_c(t);

figure(4)
subplot(3,1,1)
plot(t,y_a)
subplot(3,1,2)
plot(t,y_b)
subplot(3,1,3)
plot(t,x_c)

close all;

X = [x_a'
    y_a'
    x_b'
    y_b'
    x_c'
    y_c'];

% Code from the book:
[m,n]=size(X); % compute data size
[u,s,v]=svd((X-mean(X,2))/sqrt(n-1)); % perform the SVD
lambda=diag(s).^2; % produce diagonal variances
Y=u'*X; % produce the principal components projection

fig1 = figure(1);
sgtitle('Principle Components for Test 3')
subplot(3,1,1)
Y = Y';
plot(t, Y(:,1))
title('1st PC')
ylabel('Displacment (px)')
xlabel('Video Frame (number)')
V = axis;
V(1:2) = [0 length(t)];
axis(V)

subplot(3,1,2)
plot(t, Y(:,2))
V = axis;
V(1:2) = [0 length(t)];
axis(V)
ylabel('Displacment (px)')
xlabel('Video Frame (number)')
title('2nd PC')

subplot(3,1,3)
plot(t, Y(:,3))
V = axis;
V(1:2) = [0 length(t)];
axis(V)
ylabel('Displacment (px)')
xlabel('Video Frame (number)')
title('3rd PC')

fig2 = figure(2);
plot(1:length(lambda), lambda, '-o', 'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.75,0.75,1])
xticks(1:length(lambda))
ylabel('Variance')
xlabel('Principal Component')
title(sprintf('Magnitude of Variance Captured \n by each Principal Component \n (Test 3)'))

%% Saving

saveas(fig1,'PC3.png')
saveas(fig2,'V3.png')