%% (test 2) Noisy Case
clear variables; close all; clc;
%% Cam 1_2

load cam1_2.mat

datasize = size(vidFrames1_2); % height, width, color channels, length

x_a = zeros([datasize(4) 1]); y_a = x_a;

for j=1:datasize(4)
    gryfrm = rgb2gray(vidFrames1_2(:,:,:,j));
    frm = gryfrm>=240;
    frm(:,1:300) = 0; % remove artifact from left side
    frm(1:150,:) = 0; % remove artifact from top
    %imshow(frm, []); 
    %hold on
    [x_a(j), y_a(j)] = centroid(frm);
    %scatter(xcom, ycom, 50, 'x');
    %hold off
    %drawnow
end

clear vidFrames1_2;
'Cam1_2 done.'

%plot(y_a)

%% Cam 2_2

load cam2_2.mat;

x_b = zeros([datasize(4) 1]); y_b = x_b;

for j=1:datasize(4)
    gryfrm = rgb2gray(vidFrames2_2(:,:,:,j));
    frm = gryfrm>=245;
    frm(1:175,1:175) = 0; % remove artifact from top left side
    frm(:,end-150:end) = 0; % remove artifact from right side
    %imshow(frm, []); 
    %hold on
    [x_b(j), y_b(j)] = centroid(frm);
    %scatter(xcom, ycom, 50, 'x');
    %hold off
    drawnow
end

clear vidFrames2_2;
'Cam2_2 done.'

%plot(y_b)

%% Cam 3_2
load cam3_2.mat;

x_c = zeros([datasize(4) 1]); y_c = x_c;

for j=1:datasize(4)
    gryfrm = rgb2gray(vidFrames3_2(:,:,:,j));
    frm = gryfrm>=240;
    frm(:,1:250) = 0; % remove artifact from left side
    %imshow(frm, []); 
    %hold on
    [x_c(j), y_c(j)] = centroid(frm);
    %scatter(xcom, ycom, 50, 'x');
    %hold off
    drawnow
end

clear vidFrames3_2;
'Cam3_2 done.'

%plot(x_c)

%% (test 1) Noisy Case Analysis

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

syncpeaks = [loca(end-1) locb(end-1) locc(end-1)];

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

figure(1)
subplot(3,1,1)
Y = Y';
plot(t, (Y(:,1)-Y(1,1))/max(abs(Y(:,1))))
title('1st PC, 1st PC filtered in frequency domain')
ylabel('Normalized Displacment')
xlabel('Video Frame (number)')

Y1t = fftshift(fft(Y(:,1)));
n = length(t); % fourier modes
L = t(end) - t(1); % length of video
k=(2*pi/L)*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
figure(5)
subplot(3,1,1)
plot(ks, abs(Y1t))
xlabel('k (wave number)')
ylabel('Wave Amplitude')
title('Frequency Domain Representation of 1st PC')
axis tight
sigma = .075;
mu = .1681;
filter = 1/2 * (Gaussianfilter(ks,sigma,mu) + Gaussianfilter(ks,sigma,-mu));
subplot(3,1,2)
plot(ks,filter)
xlabel('k (wave number)')
ylabel('Wave Amplitude')
title('Mixed-Gaussian Filter (Norm 1)')
axis tight
subplot(3,1,3)
Y1tf = filter'.*Y1t;
plot(ks, abs(Y1tf))
xlabel('k (wave number)')
ylabel('Wave Amplitude')
title(sprintf('Filtered Frequency Domain Representation \n of 1st PC'))
axis tight
Y1f = ifft(ifftshift(Y1tf));
figure(1)
subplot(3,1,1)
hold on
plot(t,(Y1f-Y1f(1))/max(abs(Y1f(:))))
hold off

% V = axis;
% V(1:2) = [0 length(t)];
V = [0 length(t) -0.9 .9];
axis(V)

legend('Unfiltered','Filtered','Location','Best')

fig1 = figure(1);
sgtitle('Principle Components for Test 2')
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
title(sprintf('Magnitude of Variance Captured \n by each Principal Component \n (Test 2)'))

%% Saving

saveas(fig1,'PC2.png')
saveas(fig2,'V2.png')