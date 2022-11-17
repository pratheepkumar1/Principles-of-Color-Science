% HW 12: spectral imaging
% MJMurdoch 202010, using Roy Berns' homework and images

% load all the target images into a "stack"
img(:,:,:,1) = single(imread('slot0.tif'))./(2.^16-1);
img(:,:,:,2) = single(imread('slot1.tif'))./(2.^16-1);
img(:,:,:,3) = single(imread('slot2.tif'))./(2.^16-1);
img(:,:,:,4) = single(imread('slot3.tif'))./(2.^16-1);
img(:,:,:,5) = single(imread('slot4.tif'))./(2.^16-1);
img(:,:,:,6) = single(imread('slot5.tif'))./(2.^16-1);
img(:,:,:,7) = single(imread('slot6.tif'))./(2.^16-1);

% load all the white images into a "stack"
imgw(:,:,:,1) = single(imread('white0.tif'))./(2.^16-1);
imgw(:,:,:,2) = single(imread('white1.tif'))./(2.^16-1);
imgw(:,:,:,3) = single(imread('white2.tif'))./(2.^16-1);
imgw(:,:,:,4) = single(imread('white3.tif'))./(2.^16-1);
imgw(:,:,:,5) = single(imread('white4.tif'))./(2.^16-1);
imgw(:,:,:,6) = single(imread('white5.tif'))./(2.^16-1);
imgw(:,:,:,7) = single(imread('white6.tif'))./(2.^16-1);

%% combine RGB to monochrome by adding RGB, then flatfield

im = squeeze(sum(img,3)); % sum across the 3rd dimension (R,G,B)
imw = squeeze(sum(imgw,3));
imageff = im ./ imw;

%% show images

% uncomment if you dare: will take time and memory to display
% figure
% for i = 1:7
%    subplot(4,7,i)
%    imshow(imresize(img(:,:,:,i),1/4));
%    if i==1, ylabel('Original'); end
%    subplot(4,7,i+7)
%    imshow(imresize(imgw(:,:,:,i),1/4));
%    if i==1, ylabel('White'); end
%    subplot(4,7,i+14)
%    imshow(imresize(im(:,:,i),1/4));
%    if i==1, ylabel('Raw single channel'); end
%    subplot(4,7,i+21)
%    imshow(imresize(imageff(:,:,i),1/4)); 
%    if i==1, ylabel('FF single channel'); end
% end


%% rescale

scalar = (max(max(imageff)));
scalar = scalar.*(.5);

% scale everything together to keep color balance
image = imageff ./ max(max(max(scalar)));

% change to a 2-d variable
[nr nc nb]=size(image);
pixels=reshape(image, nr*nc, nb);

%% use mask for sampling patches (find avg camera signals for each patch)

mask = imread('mask_ccsg.tif');
mask7 = repmat(mask,1,1,7);

% image of CCSG is rotated CCW 90 degrees
% patchC will be a vector N patches by 7 channels
patchCraw = zeros(max(max(mask)),7);
for i = 1:max(max(mask))
    patchCraw(i,:) = mean(reshape(image(mask7==i),[],7));
end
% reshape to 3D array and "unrotate"
patchC3 = rot90(permute(reshape(patchCraw,[10 14 7]),[2 1 3]),3);
% convert to 2D list in CCSG order
patchC = reshape(permute(patchC3,[2 1 3]),[],7);



%% load CCSG spectral data

[num,txt] = xlsread('ccsg.xlsx');
SGref = num(:,5:end)';
SGwl = (380:10:730)';
[num,txt] = xlsread('all_1nm_data.xlsx');
cmf2 = interp1(num(:,1),num(:,6:8),SGwl);
D65 = interp1(num(:,1),num(:,3),SGwl);
SGXYZ = ((cmf2' * diag(D65) * SGref) ./ (cmf2(:,2)' * D65) )';

%% plot CCSG spectral reflectance
figure;
plot(SGwl,SGref)
xlabel('wl'); ylabel('reflectance factor');

%% Compute CCSG D65 CIELAB 
% render sRGB colors: clipped to gamut-map
SGRGB = min(max(xyz2rgb(SGXYZ),0),1);

SGLab = xyz2lab(SGXYZ);
figure;
scatter3(SGLab(:,2),SGLab(:,3),SGLab(:,1),50,SGRGB,'filled');
axis equal
xlabel('a*'); ylabel('b*'); zlabel('L*')

%% find a matrix transform from 7-channel to XYZ

% Compute XYZ matrix from CCSG camera data and known XYZ
% Use Matlab matrix-left-divide operator: \

% example M: this will provide a starting point, but is not accurate
M = [-0.5   1.6  -2.4   1.8   0.5   0.5  -0.2
     -0.4   1.0  -1.0   2.2   0.0   0.3  -0.2
     -0.3   3.7  -3.1   1.9  -0.8   0.5  -0.2 ];

% replace M with your optimized matrix

% compute Delta E (D65 white already)
SGEstXYZ = (M*patchC')';
SGEstLab = xyz2lab(SGEstXYZ);
SGDE = deltaE00(SGEstLab',SGLab')';

%% Plot difference in CIELAB

figure;
% plot actual CCSG colors
scatter3(SGLab(:,2),SGLab(:,3),SGLab(:,1),50,SGRGB,'filled');
axis equal
xlabel('a*'); ylabel('b*'); zlabel('L*')
% plot estimated CCSG colors
hold on
scatter3(SGEstLab(:,2),SGEstLab(:,3),SGEstLab(:,1),50,SGRGB);
legend('Actual','Estimate')

%% convert the whole image to XYZ and then sRGB

% go from Cam Signals to XYZ to sRGB
pixXYZ = (pixels)*M';
pixRGB = xyz2rgb(pixXYZ);
imgRGB = reshape(pixRGB,[nr nc 3]);

figure
imshow(imgRGB)

imwrite(imgRGB,'ImageOutput_sRGB.jpg');

%% spectral estimation



