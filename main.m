%   Main Function Example
%   For review
%   Input images are available at https://sites.google.com/site/photometricstereodata/

datadir = 'pmsData/ballPNG';
bitDepth = 16;
resize = 1;
gamma = 1;

data = load_datadir_re(datadir, bitDepth, resize, gamma);

image_height = size(data.imgs{1},1);
image_width = size(data.imgs{1},2);
n_scene = numel(data.imgs);

light_directions = data.s;
img = cellfun(@rgb2gray,data.imgs,'un',0);
I = cat(3,img{:});
I = reshape(I,[image_height*image_width,n_scene]);
I = I';

mask = data.mask;
mask = mask(:);

clear('data')
clear('img')

n_map = EstimateNormal(light_directions,I,0,mask);
n_map = reshape(n_map',[image_height,image_width,3]);
n_map = n_map/2 + 0.5;
imagesc(n_map)

