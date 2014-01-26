

datadir = '.\database\2obj\';
outdir = '.\database\2obj_gt\';
inputdir = '.\database\2obj_input\';
load('.\database\2obj\img_list.mat');
for i = 1:length(fls)
    filename = fls(i).name;
    openfile = [datadir filename '\' 'human_seg' '\'];
    human_names = dir([openfile '*' 'png']);
    colordir = [datadir filename '\' 'src_color' '\'];
    colorname = dir([colordir '*' 'png']);
    bwimg = [];
    for j = 1:length(human_names)
        img = imread([openfile human_names(j).name]);
        gray = rgb2gray(img);
        
        index11 = find(img(:,:,1)==255);
        index12 = find(img(:,:,2)==0);
        index13 = find(img(:,:,3)==0);
        index1 = intersect(index11,index12);
        index1 = intersect(index1,index13);
        
        index21 = find(img(:,:,3)==255);
        index22 = find(img(:,:,2)==0);
        index23 = find(img(:,:,1)==0);
        index2 = intersect(index21,index22);
        index2 = intersect(index2,index23);
        
        index = [index1;index2];
        gray(index) = 1;
        ind = 1:numel(img(:,:,1));
        gray(setdiff(ind,index)) = 0;
        bwimg(:,:,j) = gray;
    end
    bwimg = round(sum(bwimg,3)./j);
    inputimg = imread([colordir colorname(1).name]);
    imwrite(inputimg,[inputdir colorname(1).name]);
    imwrite(bwimg,[outdir colorname(1).name]);
end