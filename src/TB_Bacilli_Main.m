
close all;
clear all;
clc;

%%%% Input Image File Selection for Processing 
cd 'Image'
File = uigetfile('*.*','Select Input Image File');
mData = imread(File);      
cd ..

tic;                                                                        % Start a Stopwatch Timer
mData = imresize(mData,[480,640]);                                          % Resize an Image into fixed resolution 640x480 
        
if (ismatrix(mData))                                                        % Convert 2D matrix into 3-Channel 2D matrix if input is 2D Grayscale  
    [X,Map] = gray2ind(mData,128);
    mData = ind2rgb(X,Map);
end
rgb_inp = mData;
[Rw,Cw,Pn] = size(mData);                                      

R = mData(:,:,1);                                                           % RGB Planes Separation
G = mData(:,:,2);
B = mData(:,:,3);

%%%_Coarse Segmentation with Fixed Threshold
eR = imadjust(R,stretchlim(R));                                             % Enhancing contrast of RGB planes using Contrast Stretching
eG = imadjust(G,stretchlim(G));
eB = imadjust(B,stretchlim(B));

eI(:,:,1) = eR;                                                             % Planes Concatenation
eI(:,:,2) = eG;
eI(:,:,3) = eB;

Coarse_Out_1 = zeros(size(eR));

for i=1:1:size(eR,1)                                                        % Segmentation of interest Pixels from an image
    for j=1:1:size(eR,2)
        if ((eR(i,j)>eG(i,j)) && (eB(i,j)>=eG(i,j)) && (eG(i,j)<140) && (eR(i,j)>50) && (eR(i,j)-eG(i,j))>=30)
            Coarse_Out_1(i,j) = 1;
        end
    end
end


mform = makecform('srgb2lab');
Lab = applycform(mData,mform);

L = Lab(:,:,1);                                                             % Extracting Luminance component from an Image

Win_size = 31;
k_sau = 0.3;

Idata = double(L);
Mxy = medfilt2(Idata,[Win_size,Win_size]); 

Var_i = ((Idata - mean2(Idata)).^2);
Std_i = sqrt(Var_i);

Sauv_temp = (Std_i./max(Std_i(:))) - 1;
T_sauv = Mxy.*(1+(k_sau.*Sauv_temp));                                       

Bxy_sav = Idata > T_sauv;
Bxy_sav = imcomplement(Bxy_sav);

Bsav_1 = bwareaopen(Bxy_sav,100);                                           % Remove large unwanted objects from background
Bsav_2 = Bxy_sav - Bsav_1;

Bsav_3 = bwareaopen(Bsav_2,10,8);
Dpixs = bitxor(Coarse_Out_1,Bsav_3);
Diff_1 = logical(abs(Bsav_3 - Dpixs));
Seg = bitor(Coarse_Out_1,Diff_1);                                           % Merge Coarse segmentation with Fine segmentation output

% Retain bacilli pixels and Filter Non bacilli pixels from segmented Image 
[S_label,S_cnt] = bwlabel(Seg,8);                                           % 8-pixel connected component analysis and find shape descriptors  
S_props = regionprops(S_label,'All');

B_area = zeros(1,length(S_props));
B_Peri = zeros(1,length(S_props));
B_Hperi = zeros(1,length(S_props));

for Bi=1:1:length(S_props)
    B_area(1,Bi) = S_props(Bi).Area;                % Area
    B_Peri(1,Bi) = S_props(Bi).Perimeter;           % Perimeter
    B_cImage = S_props(Bi).ConvexImage;
    H_prop = regionprops(B_cImage,'Perimeter');
    B_Hperi(1,Bi) = H_prop.Perimeter;               % Convex Hull Image Perimeter
end

Roughness = B_Peri./B_Hperi;                        % Roughness
Circ = ((4*pi).*B_area)./B_Peri.^2;                 % Circularity
Major = [S_props.MajorAxisLength];                  % Major Axis Length
        
%%%Remove Non Objects from Background and Keep Objects
A_Th1 = 200; A_Th2 = 5; R_Th = 1.5;
Max_Th1 = 4; Max_Th2 = 28;
S_Out = S_label;

Bac_Obj_A = (B_area > A_Th1) | (B_area < A_Th2);            % Area based Filter
A_ind = find(Bac_Obj_A);

for ii=1:1:length(A_ind)
    S_Out(S_Out==A_ind(ii)) = 0;
end

Bac_Obj_Major = (Major < Max_Th1) | (Major > Max_Th2);      % Major axis length based Filter
M_ind = find(Bac_Obj_Major);

for ii=1:1:length(M_ind)
    S_Out(S_Out==M_ind(ii)) = 0;
end

Bac_Obj_R = (Roughness >= R_Th);                            % Roughness based Filter
R_ind = find(Bac_Obj_R);

for ii=1:1:length(R_ind)
    S_Out(S_Out==R_ind(ii)) = 0;
end

[S_Out_label,S1_cnt] = bwlabel(S_Out,8);                                    % Find Number of Objects and Bounding Box Vector 
S_props_1 = regionprops(S_Out_label,'BoundingBox');
BBox = [S_props_1.BoundingBox];
BBox_1 = (reshape(BBox,[4,length(BBox)/4]))';
E_time = toc;                                                               

rgb_inp = insertShape(rgb_inp,'Rectangle',BBox_1 ,'Color','Red');           % Draw Rectangle shape box on each segmented object
rgb_inp = insertText(rgb_inp,[Cw-180,Rw-30],['Count(Approx) - ',num2str(S1_cnt)],'FontSize',10,'TextColor','Red','BoxColor','White');  % Annotate Masured Bacilli Count

%%% Display Measured TB bacilli Count-Execution time-Processed Images
fprintf(1,'%s %s\n','Processed File -',File);
fprintf(1,'%s\n','---------------------------------');
fprintf(1,'%s - %d\n','Count(Approx)',S1_cnt);
fprintf(1,'%s\n','---------------------------------');
fprintf(1,'%s : %f\n','Elapsed Time(sec)',E_time);

figure('Name','Input Image','MenuBar','None');
imshow(mData);

figure('Name','Segmented Object','MenuBar','None');
imshow(Seg);

figure('Name','Postprocessed Segmented Image','MenuBar','None');
imshow(S_Out);

figure('Name','Detected Object','MenuBar','None');
imshow(rgb_inp);

