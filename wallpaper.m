% Title: Wallpaper GUI, 2023 Summer Research Project (Part II)
% By Paisley (Jingxuan) Hou
% Date: 07/21/2023
% ------Description------
% This program creates a GUI that displays wallpaper patterns along with
% its original image based on a user-defined symmetry group. The program
% incorporates options for all 15 possible symmetry groups for the wallpapers.
% The user can also specify the number of terms they desire for the wave function,
% and any RGB image from which to generate the wallpaper based on the domain-coloring
% algorithm provided in the program.
% The user interface contains sliders for users to play around with 
% different coefficients and frequencies for the wave function and
% generate different wallpaper patterns from the given image.
% Mathematical theories behind the program referred to the book Creating Symmetry, by Professor Farris from Santa Clara University.


%INPUTs

%User-defined inputs
sym = input("Symmetry group: ", "s");
nTerms = input("Number of terms: ");
img_path = input("Image file: ", "s");

img = imread(img_path);
%%
%Define dual lattices
if sym == "p2"
    a = 1/2;
    b = sqrt(3)/2;
    k1 = 1;
    k2 = complex(a, b);
elseif sym == "p3" | sym == "p31m" | sym == "p3m1" | sym == "p6" | sym == "p6m"
    k1 = 1;
    k2 = complex(-1/2, sqrt(3)/2);
elseif sym == "p4" | sym == "p4m" | sym == "p4g" | sym == "p4g/cmm"
    k1 = 1;
    k2 = i;
elseif sym == "cm" | sym == "cmm"
    b = sqrt(3)/2;
    k1 = complex(1/2, b);
    k2 = complex(1/2, -b);
elseif sym == "pm" | sym == "pmm" |sym == "pmg" | sym == "pg" | sym == "pgg" 
    L = 1.5;
    k1 = 1;
    k2 = complex(0,L);
else
    disp("Invalid symmetry group.");
    exit
end
%%
%Set up graphical interface

%Set up figure
fig = uifigure("WindowStyle","normal","Name","Wallpaper("+sym+")","Position",[100 100 500 500]);

%Customize figure layout
g = uigridlayout(fig,[nTerms*2+1,8], "BackgroundColor","#ffc0cb");
g.RowHeight(1) = {20};
rowHeight = {480/(2*nTerms)};
g.RowHeight(2:nTerms*2+1) = repmat(rowHeight,1,nTerms*2);
g.ColumnWidth = {100,100,100,100,'3x','3x','3x','3x'};

%Set up axis for images

%Rosette Image
ax_wave = uiaxes(g);
ax_wave.Layout.Row = [2 nTerms*2+1];
ax_wave.Layout.Column = [5 6];

%Original Image
ax_orig = uiaxes(g);
ax_orig.Layout.Row = [2 nTerms*2+1];
ax_orig.Layout.Column = [7 8];
%%
%Set buttons
showBtn = uibutton(g, "push", "Text", "Show Wallpaper", "ButtonPushedFcn", @(src,event) showWallpaper(ax_wave));
showBtn.Layout.Row = 1;
showBtn.Layout.Column = 6;
%%
%Set labels

%coef
coefLabel = uilabel(g,"Text","coef.(r,θ)");
coefLabel.HorizontalAlignment = 'center';
coefLabel.Layout.Row = 1;
coefLabel.Layout.Column = 1;
coefLabel.FontName = "TimesNewRoman";
coefLabel.FontSize = 10;

%freq
freqLabel = uilabel(g,"Text","freq.(n&m)");
freqLabel.HorizontalAlignment = 'center';
freqLabel.Layout.Row = 1;
freqLabel.Layout.Column = 3;
freqLabel.FontName = "TimesNewRoman";
freqLabel.FontSize = 10;

%domain

%zmin
zminLabel = uilabel(g, "Text", "zMin");
zminLabel.HorizontalAlignment = 'left';
zminLabel.Layout.Row = 1;
zminLabel.Layout.Column = 4;
zminLabel.FontName = "TimesNewRoman";
zminLabel.FontSize = 10;

%zmax
zmaxLabel = uilabel(g, "Text", "zMax");
zmaxLabel.HorizontalAlignment = 'left';
zmaxLabel.Layout.Row = 3;
zmaxLabel.Layout.Column = 4;
zmaxLabel.FontName = "TimesNewRoman";
zmaxLabel.FontSize = 10;

%zmin
znumLabel = uilabel(g, "Text", "zNum");
znumLabel.HorizontalAlignment = 'left';
znumLabel.Layout.Row = 5;
znumLabel.Layout.Column = 4;
znumLabel.FontName = "TimesNewRoman";
znumLabel.FontSize = 10;

%wallpaper
waveLabel = uilabel(g,"Text","Wallpaper", "FontWeight","bold");
waveLabel.HorizontalAlignment = 'center';
waveLabel.Layout.Row = 1;
waveLabel.Layout.Column = 5;
waveLabel.FontName = "TimesNewRoman";
waveLabel.FontSize = 10;

%original image
imgLabel = uilabel(g,"Text","Original Image", "FontWeight","bold");
imgLabel.HorizontalAlignment = 'center';
imgLabel.Layout.Row = 1;
imgLabel.Layout.Column = [7 8];
imgLabel.FontName = "TimesNewRoman";
imgLabel.FontSize = 10;
%%
%Inialize window display
initializeFigure(nTerms, sym, k1, k2, ax_wave, ax_orig, img);
create_controllers(g, nTerms, sym, k1,k2, img, ax_wave);
%%
function initializeFigure(nTerms, sym, k1, k2, ax_wave, ax_orig, img)
%Initialize the figure

%Initiliaze slider-controlled inputs

global Magnitude Theta Coef; %coefficients
global nFreq mFreq; %frequencies
global zMin zMax zNum; %domain

%Initialize coefficients
Magnitude = ones(nTerms,1);
Theta = zeros(nTerms,1);
Coef = complex(Magnitude.*cos(Theta), Magnitude.*sin(Theta));

%Initialize frequencies
mFreq = ones(nTerms,1);
nFreq = zeros(nTerms,1);

%Initialize domain
zMin = -2;
zMax = 2;
zNum = 100;

%Set up function
f_val = generate_func(sym,k1,k2);

%display rosette image
domain_coloring(img, f_val, ax_wave);

%display orignal image
imshow(img, "Parent", ax_orig);

end


function f_val = generate_func(sym,k1,k2)
%Generate function values
global Coef
global nFreq mFreq
global zMin zMax zNum

[x, y] = meshgrid(linspace(zMin, zMax, zNum));
z = complex(x,y);
[X, Y] = generate_XY(k1,k2,z);

f_val = zeros(zNum,zNum);
nTerms = size(Coef,1);

if sym == "p3"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/3 * ( exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(mFreq(k) * X - (nFreq(k) + mFreq(k)) * Y)) + ...
            exp(2*pi*i*(-(nFreq(k)+mFreq(k)) * X + nFreq(k) * Y)));
    end

elseif sym == "p31m"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/3 * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(mFreq(k) * X - (nFreq(k) + mFreq(k)) * Y)) + ...
            exp(2*pi*i*(-(nFreq(k)+mFreq(k)) * X + nFreq(k) * Y)) + ...
            exp(2*pi*i*(mFreq(k) * X + nFreq(k) * Y)) + ...
            exp(2*pi*i*(nFreq(k) * X - (mFreq(k) + nFreq(k)) * Y)) + ...
            exp(2*pi*i*(-(mFreq(k)+nFreq(k)) * X + mFreq(k) * Y)));
    end

elseif sym == "p3m1"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/3 * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(mFreq(k) * X - (nFreq(k) + mFreq(k)) * Y)) + ...
            exp(2*pi*i*(-(nFreq(k)+mFreq(k)) * X + nFreq(k) * Y)) + ...
            exp(2*pi*i*(-mFreq(k) * X - nFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X + (mFreq(k) + nFreq(k)) * Y)) + ...
            exp(2*pi*i*((mFreq(k)+nFreq(k)) * X - mFreq(k) * Y)));
    end

elseif sym == "p6"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/3 * ( exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(mFreq(k) * X - (nFreq(k) + mFreq(k)) * Y)) + ...
            exp(2*pi*i*(-(nFreq(k)+mFreq(k)) * X + nFreq(k) * Y)) + ...
            exp(2*pi*i*((-nFreq(k))*X + (-mFreq(k))*Y)) + ...
            exp(2*pi*i*((-mFreq(k))*X + (nFreq(k)+mFreq(k))*Y)) + ...
            exp(2*pi*i*((nFreq(k)+mFreq(k))*X + (-nFreq(k))*Y)));
    end

elseif sym == "p6m"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/3 * ( exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(mFreq(k) * X - (nFreq(k) + mFreq(k)) * Y)) + ...
            exp(2*pi*i*(-(nFreq(k)+mFreq(k)) * X + nFreq(k) * Y)) + ...
            exp(2*pi*i*((mFreq(k))*X + (nFreq(k))*Y)) + ...
            exp(2*pi*i*((-(nFreq(k)+mFreq(k)))*X + (mFreq(k))*Y)) + ...
            exp(2*pi*i*((nFreq(k))*X + (-(nFreq(k)+mFreq(k)))*Y)) + ...
            exp(2*pi*i*((-nFreq(k))*X + (-mFreq(k))*Y)) + ...
            exp(2*pi*i*((-mFreq(k))*X + (nFreq(k)+mFreq(k))*Y)) + ...
            exp(2*pi*i*((nFreq(k)+mFreq(k))*X + (-nFreq(k))*Y)) + ...
            exp(2*pi*i*((-mFreq(k))*X + (-nFreq(k))*Y)) + ...
            exp(2*pi*i*(((nFreq(k)+mFreq(k)))*X - (mFreq(k))*Y)) + ...
            exp(2*pi*i*((-nFreq(k))*X + ((nFreq(k)+mFreq(k)))*Y)));
    end

elseif sym == "p4"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/4 * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(mFreq(k) * X - nFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X - mFreq(k) * Y)) + ...
            exp(2*pi*i*(-mFreq(k) * X + nFreq(k) * Y)) ...
            );
    end

elseif sym == "p4m"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/4 * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(mFreq(k) * X - nFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X - mFreq(k) * Y)) + ...
            exp(2*pi*i*(-mFreq(k) * X + nFreq(k) * Y)) + ...
            exp(2*pi*i*(mFreq(k) * X + nFreq(k) * Y)) + ...
            exp(2*pi*i*(nFreq(k) * X - mFreq(k) * Y)) + ...
            exp(2*pi*i*(-mFreq(k) * X - nFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X + mFreq(k) * Y)) ...
            );
    end

elseif sym == "p4g"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/4 * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(mFreq(k) * X - nFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X - mFreq(k) * Y)) + ...
            exp(2*pi*i*(-mFreq(k) * X + nFreq(k) * Y)) + ...
            (-1).^(nFreq(k) + mFreq(k)) * ( ...
            exp(2*pi*i*(mFreq(k) * X + nFreq(k) * Y)) + ...
            exp(2*pi*i*(nFreq(k) * X - mFreq(k) * Y)) + ...
            exp(2*pi*i*(-mFreq(k) * X - nFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X + mFreq(k) * Y)) ...
            ));
    end

elseif sym == "cm"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/2 * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(mFreq(k) * X + nFreq(k) * Y)) ...
            );
    end

elseif sym == "cmm"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/2 * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(mFreq(k) * X + nFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X - mFreq(k) * Y)) + ...
            exp(2*pi*i*(-mFreq(k) * X - nFreq(k) * Y)) ...
            );
    end

elseif sym == "pm"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X + mFreq(k) * Y)) ...
            );
    end

elseif sym == "pmm"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/2 * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X - mFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(nFreq(k) * X - mFreq(k) * Y)) ...
            );
    end

elseif sym == "pmg"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/2 * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X - mFreq(k) * Y)) + ...
            (-1).^nFreq(k) * (exp(2*pi*i*(nFreq(k) * X - mFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X + mFreq(k) * Y))) ...
            );
    end

elseif sym == "pg"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            (-1)^mFreq(k) * exp(2*pi*i*(-nFreq(k) * X + mFreq(k) * Y)) ...
            );
    end

elseif sym == "pgg"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * 1/2 * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X - mFreq(k) * Y)) + ...
            (-1).^(nFreq(k)+mFreq(k)) * (exp(2*pi*i*(nFreq(k) * X - mFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X + mFreq(k) * Y))) ...
            );
    end

elseif sym == "p4g/cmm"
    for k = 1:nTerms
        f_val = f_val + Coef(k) * ( ...
            exp(2*pi*i*(nFreq(k) * X + mFreq(k) * Y)) - ...
            exp(2*pi*i*(mFreq(k) * X - nFreq(k) * Y)) + ...
            exp(2*pi*i*(-nFreq(k) * X - mFreq(k) * Y)) - ...
            exp(2*pi*i*(-mFreq(k) * X + nFreq(k) * Y)) + ...
            (-1).^(nFreq(k) + mFreq(k)) .* ( ...
            exp(2*pi*i*(mFreq(k) * X + nFreq(k) * Y)) - ...
            exp(2*pi*i*(nFreq(k) * X - mFreq(k) * Y)) + ...
            exp(2*pi*i*(-mFreq(k) * X - nFreq(k) * Y)) - ...
            exp(2*pi*i*(-nFreq(k) * X + mFreq(k) * Y)) ...
            ));
    end

end

end


function [X,Y] = generate_XY(k1,k2,z)
v1 = (i*k2)/imag(k1*conj(k2));
v2 = (i*k1)/imag(conj(k1)*k2);
X = real(z*conj(v1));
Y = real(z*conj(v2));
end


function domain_coloring(img, f_val, ax_wave)
%Creates the domain-colored image from the given complex-valued function
%and draws it on the given axis ax_wave.

%Share global variables
global Coef nFreq mFreq

%Get dimensions of the reference image
[nRows, nCols] = size(img,1:2);
n = min(nRows, nCols);

%Get dimensions of the function space
zNum = size(f_val, 1);

%Calculate polar form of function space
Radius = abs(f_val);
Angle = angle(f_val);

%Scale radius to be within the range of the image size
scaledRadius = scale(Radius,0,n/2);

%calculate x and y coords
X_coords = scaledRadius .* cos(Angle);
Y_coords = scaledRadius .* sin(Angle);

%Scale x and y coords to match the pixel positions on the image
%Uses the central square section of the image
if nRows < nCols
    X_coords_scaled = round(scale(X_coords, floor(nCols/2 - n/2 + 0.5), ceil(nCols/2 + n/2 + 0.5)));
    Y_coords_scaled = n + 1 - round(scale(Y_coords, 1, n));
else
    X_coords_scaled = round(scale(X_coords, 1, n));
    Y_coords_scaled = round(scale(-Y_coords, floor(nRows/2 - n/2 + 0.5), ceil(nRows/2 + n/2 + 0.5)));
end


%Fill in new RGB matrix
newImg = zeros(zNum,zNum,3,"uint8");
for i = 1:zNum
    for j = 1:zNum
        x = X_coords_scaled(i,j);
        y = Y_coords_scaled(i,j);
        newImg(i,j,:) = img(y,x,:);
    end
end

%display new image
imshow(newImg, "Parent", ax_wave);

%display corresponding variable values
values = "";

nTerms = size(Coef, 1);

for k = 1:nTerms
    %coef
    values = values + Coef(k) + ", ";
    %freq
    values = values + nFreq(k) + ", ";
    values = values + mFreq(k) + "|";
end

title(ax_wave, values, FontSize=7);

end


function create_controllers(g, nTerms, sym, k1,k2, img, ax_wave)
for i = 1:nTerms
    %coef sliders

    %magnitude
    sld = uislider(g,"Value",1, "Limits", [0 50], "ValueChangingFcn", ...
        @(src,event)updateMagnitude(src,event,i, nTerms, sym, k1,k2, img, ax_wave));
    sld.Layout.Row = 2*i;
    sld.Layout.Column = 1;

    %angle
    gauge = uigauge(g, "BackgroundColor","#c0cbff","ScaleDirection","counterclockwise", ...
        "Limits", [0, 2*pi], "MajorTicks", [0:pi/2:2*pi], "MajorTickLabels", {'0','π/2','π','3π/2','2π'}, ...
        "MinorTicks", [], "Value", 0);
    gauge.Layout.Row = 2*i+1;
    gauge.Layout.Column = 2;

    sld = uislider(g,"Value",0,"Limits", [0 2*pi], "ValueChangingFcn",...
        @(src,event)updateTheta(src,event,i, gauge, nTerms, sym, k1,k2, img, ax_wave));
    sld.Layout.Row = 2*i+1;
    sld.Layout.Column = 1;

    %freq sliders

    %n
    sld = uislider(g,"Value",0, "Limits", [-20 20], "ValueChangingFcn", ...
        @(src,event)update_n(src, event,i,nTerms, sym, k1,k2, img, ax_wave));
    sld.Layout.Row = 2*i;
    sld.Layout.Column = 3;

    %m
    sld = uislider(g,"Value",1, "Limits", [-20 20], "ValueChangingFcn", ...
        @(src,event)update_m(src, event,i,nTerms, sym, k1,k2, img, ax_wave));
    sld.Layout.Row = 2*i+1;
    sld.Layout.Column = 3;

end

%domain sliders

%zMin
sld = uislider(g, "Value", -2, "Limits", [-10 -1], "ValueChangingFcn", ...
    @(src, event)update_zmin(src, event, nTerms, sym, k1,k2, img, ax_wave));
sld.Layout.Row = 2;
sld.Layout.Column = 4;

%zMax
sld = uislider(g, "Value", 2, "Limits", [1 10], "ValueChangingFcn", ...
    @(src, event)update_zmax(src, event, nTerms, sym, k1,k2,img, ax_wave));
sld.Layout.Row = 4;
sld.Layout.Column = 4;

%zNum
sld = uislider(g, "Value", 100, "Limits", [100 3000], "ValueChangingFcn", ...
    @(src, event)update_znum(src, event, nTerms, sym, k1,k2,img, ax_wave));
sld.Layout.Row = 6;
sld.Layout.Column = 4;
end


function updateMagnitude(src,event,i, nTerms, sym, k1,k2, img, ax_wave)
global Magnitude Theta Coef
Magnitude(i) = event.Value;
Coef = complex(Magnitude.*cos(Theta), Magnitude.*sin(Theta));
f_val = generate_func(sym,k1,k2);
domain_coloring(img, f_val, ax_wave);
end


function updateTheta(src,event,i, gauge, nTerms, sym, k1,k2,img, ax_wave)
global Magnitude Theta Coef
Theta(i) = event.Value;
gauge.Value = event.Value;
Coef = complex(Magnitude.*cos(Theta), Magnitude.*sin(Theta));
f_val = generate_func(sym,k1,k2);
domain_coloring(img, f_val, ax_wave);
end


function update_n(src, event, i, nTerms, sym, k1,k2,img, ax_wave)
global nFreq
nFreq(i) = round(event.Value);
f_val = generate_func(sym,k1,k2);
domain_coloring(img, f_val, ax_wave);
end


function update_m(src, event, i, nTerms, sym, k1,k2, img, ax_wave)
global mFreq
mFreq(i) = round(event.Value);
f_val = generate_func(sym,k1,k2);
domain_coloring(img, f_val, ax_wave);
end


function update_zmin(src, event, nTerms, sym, k1,k2, img, ax_wave)
global zMin
zMin = event.Value;
f_val = generate_func(sym,k1,k2);
domain_coloring(img, f_val, ax_wave);
end


function update_zmax(src, event, nTerms, sym, k1,k2, img, ax_wave)
global zMax
zMax = event.Value;
f_val = generate_func(sym,k1,k2);
domain_coloring(img, f_val, ax_wave);
end


function update_znum(src, event, nTerms, sym, k1,k2, img, ax_wave)
global zNum
zNum = round(event.Value);
f_val = generate_func(sym,k1,k2);
domain_coloring(img, f_val, ax_wave);
end


function showWallpaper(ax_wave)
wallpaper = getimage(ax_wave);
figure, imshow(wallpaper);
end


function Y = scale(X, new_min, new_max)
Xmin = min(X,[],"all");
Xmax = max(X,[],"all");
if Xmax==Xmin %X values are the same
    if abs(Xmin-new_min) < abs(Xmin-new_max) %X closer to new_min
        Y = X-X + new_min;
    else
        Y = X-X + new_max;
    end
else
    Y = new_min + ((X-Xmin)/(Xmax-Xmin))*(new_max-new_min);
end

end
