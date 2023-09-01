% Title: Bessel Wave Surface GUI, 2023 Summer Research Project (Part III)
% By Paisley (Jingxuan) Hou------
% Date: 08/12/2023------
% Description:
% This program file creates a GUI window that displays the wave
% surface from bessel functions as a superposition of a user-defined number
% of terms. Each term of the wave surface is defined as follows:
% f(r,) = J_n(_{n,m}*r) * (A*(cos(n(+)) + sin(n(+)))), where
% J_n(_{n,m}*r) is the bessel function of the nth order with the mth zero.
%
% Domain-coloring is applied to make the rendered image has a more
% interesting surface. The source image for domain-coloring is given by the user. 
% The colormap generated from the source image is shown on the rightmost
% section of the interface, which can be displayed separately by clicking
% the "Show Colormap" button. The colormap is in polar coordinates, and is
% generated based on the values of the following function over the (R,)
% domain: a_{n,m} * e^(in) * J_n(_{n,m}*r), where a_{n,m} can be any
% number from the complex space.
%
% The GUI contains interactive components like sliders and spinners for users to 
% play with different parameter values for the wave function.
% For efficiency concerns, the size of the domain has been initialized to be
% 100. However, if the user wants a better resolution of the image, they can adjust
% the size using the given slider in the interface.  
%
% The user may click on the "Show Wave Surface" button to gain a clearer
% view of the wave surface in a separate window. 
%
% Another feature of the GUI is that it allows users to observe the
% animated wave surface by clicking on the "animate" button above the
% figure. The animation will be displayed in a pop-up window, and will be
% saved to a .avi file titled "vibrating_drumhead.avi" in the same
% directory as the where the current program is located.
%


%User-defined variables
nTerms = input("Number of terms: ");
imgPath = input("Source image file for texture mapping: ","s");
img = imread(imgPath);
%%
%Set up GUI window
fig = uifigure("WindowStyle","normal","Name","BesselWaveSurface","Position",[200 100 700 500]);

%Define layout
g = uigridlayout(fig, [2*nTerms+1, 9], "BackgroundColor", "#fff2cc");
g.ColumnWidth = {100, 100, 100, 100,100,'1x', '1x', '1x', '1x'};
g.RowHeight(1) = {20};
h = {480/(2*nTerms)};
g.RowHeight(2:2*nTerms+1) = repmat(h,1,2*nTerms);
%%
%Define axis for the surface plot
ax = uiaxes(g);
ax.Layout.Row = [2 2*nTerms+1];
ax.Layout.Column = [6 7];
ax.Color = "#FFE4E1";
grid(ax, "on");

%Define axis for colormap from domain-coloring
ax_2dmap = uiaxes(g);
ax_2dmap.Layout.Row = [2 2*nTerms+1];
ax_2dmap.Layout.Column = [8 9];
%%
%Set buttons

%Set button to animate the current wave surface
animateBtn = uibutton(g, "push", "Text", "Animate", "ButtonPushedFcn", @(src, event) animateWave(ax_2dmap));
animateBtn.Layout.Row = 1;
animateBtn.Layout.Column = 7;

%Set button to display surface plot in a separate window
surfShowBtn = uibutton(g, "Push", "Text", "Show Wave Surface", "ButtonPushedFcn", @(src,event) SurfShowButtonPushed(nTerms, ax_2dmap));
surfShowBtn.Layout.Row = 1;
surfShowBtn.Layout.Column = 6;

%Set button to display colormap in a separate window
cmapcmapShowBtn = uibutton(g, "push", "Text", "Show Colormap", "ButtonPushedFcn", @(src,event) CMapShowButtonPushed(ax_2dmap));
cmapcmapShowBtn.Layout.Row = 1;
cmapcmapShowBtn.Layout.Column = [8 9];
%%
%Set labels
nLabel = uilabel(g, "Text", "n", "FontSize", 15, "FontWeight", "bold");
nLabel.Layout.Row = 1;
nLabel.Layout.Column = 1;

mLabel = uilabel(g, "Text", "m", "FontSize", 15, "FontWeight", "bold");
mLabel.Layout.Row = 1;
mLabel.Layout.Column = 2;

ampLabel = uilabel(g, "Text", "Amp.(mag.&angle)", "FontSize", 3, "FontWeight", "bold");
ampLabel.Layout.Row = 1;
ampLabel.Layout.Column = 3;

numLabel = uilabel(g, "Text", "Domain Size", "FontSize", 15, "FontWeight", "bold");
numLabel.Layout.Row = 1;
numLabel.Layout.Column = 5;
%%
%Initialize GUI
initFig(g, nTerms, ax, ax_2dmap, img);
%%
function initFig(g, nTerms, ax, ax_2dmap, img)
global ns ms mags angles amps Num R Theta

%Initialize parameters
ns = zeros(nTerms,1);
ms = ones(nTerms, 1);
mags = zeros(nTerms,1);
angles = zeros(nTerms, 1);
amps = complex(mags.*cos(angles), mags.*sin(angles));

%Initialize domain size to be 100 and define domain
Num = 100;
r = 0:1/Num:1;
theta = 0:2*pi/Num:2*pi;
[R, Theta] = meshgrid(r, theta);

%Evaluate wave function values
z = zEval();

%Generate colormap for the current wave function
colormap = BesselDomainColoring(nTerms, img, ax_2dmap);

%Create surface with the evaluated function values and the generated colormap
plotSurf(ax,colormap,z,nTerms);

%Create controllers for the GUI
createControllers(g,ax,ax_2dmap, img, nTerms);

end


function z = zEval()
global ns ms mags angles Num R Theta

%Initialize function space
z = zeros(Num+1, Num+1);

%Calculate wave function values
nTerms = size(ns,1);
for i = 1:nTerms
    n = ns(i);
    m = ms(i);
    A = mags(i);
    phi = angles(i);
    L = besselzero(n, 100, 1);
    z = z + besselj(n, L(m)*R) .* (A*(cos(n.*(Theta+phi))+sin(n.*(Theta+phi))));
end

end


function plotSurf(ax,colormap,z,nTerms)
global ns ms mags angles R Theta

%Compute cartesian coordinates
x = R.*cos(Theta);
y = R.*sin(Theta);

%Find z limits
zMax = max(z, [], 'all');
zMin = -zMax;

%Plot surface
surf(ax,x,y,z,colormap, "EdgeColor", "none", "FaceColor", "interp");
zlim([zMin-5, zMax+5]);
axis manual

%Show wave function
eq = "f(r,) = ";
for i = 1:nTerms
    newTerm = "";
    n = ns(i);
    m = ms(i);
    A = mags(i);
    phi = angles(i);
    newTerm = "J_" + n + "(_" + n + ",_"+ m+"r)"+"("+A+"cos(+"+phi+")"+"+sin(+"+phi+")"+")";
    eq = eq + newTerm;
    if i < nTerms
        eq = eq + "+";
    end
end
title(ax, eq,"FontWeight","bold", "FontAngle",'italic','FontSize',9);

end


function createControllers(g,ax,ax_2dmap, img, nTerms)
global ns ms mags angles Num

%Create controllers for the set of parameters for each term
for i = 1:nTerms
    mag = mags(i);
    angle = angles(i);
    n = ns(i);
    m = ms(i);

    %Create spinners for n and m
    
    %n
    spinner = uispinner(g, "Value", n, "Limits", [0 30], "RoundFractionalValues","on", ...
        "BackgroundColor", "#f0b1a4", ...
        "ValueChangingFcn", @(spn, event)nChanged(spn, event, i, ax, ax_2dmap, img, nTerms));
    spinner.Layout.Row = 2*i;
    spinner.Layout.Column = 1;
    
    %m
    spinner = uispinner(g, "Value", m, "Limits", [1 12], "RoundFractionalValues","on", ...
        "BackgroundColor", "#FE7F8A", ...
        "ValueChangingFcn", @(spn, event)mChanged(spn, event, i, ax, ax_2dmap, img, nTerms));
    spinner.Layout.Row = 2*i;
    spinner.Layout.Column = 2;

    %Create sliders for mag. and angle for the amplitude
    
    %mag
    sld = uislider(g, "Value", mag, "Limits", [0 50], ...
        "ValueChangingFcn", @(src,event)magChanged(src,event, i, ax, ax_2dmap, img, nTerms));
    sld.Layout.Row = 2*i;
    sld.Layout.Column = 3;

    %angle

    %Create a gauge that helps visualize the postion of the angle
    gauge = uigauge(g, "BackgroundColor","#c0cbff","ScaleDirection","counterclockwise", ...
    "Limits", [0, 2*pi], "MajorTicks", 0:pi/2:2*pi, "MajorTickLabels", {'0','/2','','3/2','2'}, ...
    "MinorTicks", [], "Value", angle);
    gauge.Layout.Row = 2*i+1;
    gauge.Layout.Column = 4;

    %Slider for the angle
    sld = uislider(g, "Value", angle, "Limits", [0 2*pi], ...
        "ValueChangingFcn", @(src,event)angleChanged(src,event, i, ax, ax_2dmap, img, nTerms, gauge));
    sld.Layout.Row = 2*i+1;
    sld.Layout.Column = 3;

end

%Create sliders for domain size
sld = uislider(g, "Value", Num, "Limits", [100 3000], ...
    "ValueChangedFcn", @(src, event)numChanged(src, event, ax, ax_2dmap, img, nTerms));
sld.Layout.Row = 3;
sld.Layout.Column = 5;

end


function nChanged(spn, event, i, ax, ax_2dmap, img, nTerms)
global ns

%Get new n value
ns(i) = event.Value;

%Evaluate function for the new n value
z = zEval();

%Generate colormap for the new n value
colormap = BesselDomainColoring(nTerms, img, ax_2dmap);

%Plot surface for the new n value
plotSurf(ax, colormap, z, nTerms);

end


function mChanged(spn, event, i, ax, ax_2dmap, img, nTerms)
global ms

%Get new m value
ms(i) = event.Value;

%Evaluate function for the new m value
z = zEval();

%Generate colormap for the new m value
colormap = BesselDomainColoring(nTerms, img, ax_2dmap);

%Plot surface for the new m value
plotSurf(ax, colormap, z, nTerms);

end


function magChanged(spn, event, i, ax, ax_2dmap, img, nTerms)
global mags angles amps

%Get new mag value
mags(i) = event.Value;

%Change the corresponding amplitude
amps(i) = complex(mags(i)*cos(angles(i)), mags(i)*sin(angles(i)));

%Evaluate function for the new mag value
z = zEval();

%Generate colormap for the new mag value
colormap = BesselDomainColoring(nTerms, img, ax_2dmap);

%Plot surface for the new mag value
plotSurf(ax, colormap, z, nTerms);

end


function angleChanged(spn, event, i, ax, ax_2dmap, img, nTerms, gauge)
global angles mags amps

%Get new angle value
angles(i) = event.Value;

%Change the corresponding amplitude
amps(i) = complex(mags(i)*cos(angles(i)), mags(i)*sin(angles(i)));

%Update gauge to match the new angle value
gauge.Value = event.Value;

%Evaluate function for the new angle value
z = zEval();

%Generate colormap for the new angle value
colormap = BesselDomainColoring(nTerms, img, ax_2dmap);

%Plot surface for the new angle value
plotSurf(ax, colormap, z, nTerms);

end


function numChanged(src, event, ax, ax_2dmap, img, nTerms)
global Num R Theta

%Redefine domain with the new domain size
Num = round(event.Value);
r = 0:1/Num:1;
theta = 0:2*pi/Num:2*pi;
[R, Theta] = meshgrid(r, theta);

%Evaluate function on the new domain
z = zEval();

%Generate colormap with the number of pixels on each side matching the new
%domain size
colormap = BesselDomainColoring(nTerms, img, ax_2dmap);

%Plot surface on the new domain
plotSurf(ax, colormap, z, nTerms);

end


function CMap = BesselDomainColoring(nTerms, img, ax_2dmap)
global ns ms amps R Theta

%Evaluate function values
Num = size(R, 1);
z = zeros(Num, Num);
for k = 1:nTerms
    n = ns(k);
    m = ms(k);
    amp = amps(k);
    L = besselzero(n, 20, 1);
    J = besselj(n, L(m).*R);
    z = z + amp.*J.*complex(cos(n*Theta), sin(n*Theta));
end
RealPart = real(z);
ImagPart = imag(z);

%Rescale function values to match the (x,y) coordinates on
%the original image
[nRows, nCols] = size(img, 1:2);
RealPartScaled = round(scale(RealPart, 1, nCols));
ImagPartScaled = nRows + 1 - round(scale(ImagPart, 1, nRows));

%Initialize CMap matrix
CMap = zeros(Num,Num,3,"uint8");

%Domain Coloring
for i = 1:Num
    for j = 1:Num
        %Get real, imag. parts of the corresponding function value.
        %Each (real, imag.) pair maps to a pixel on the given image.
        x = RealPartScaled(i,j);
        y = ImagPartScaled(i,j);

        %Fill in the corresponding pixel of colormap
        CMap(i,j,:) = img(y,x,:);
    end
end


%Display colormap
CMapRot = imrotate(CMap, 270);
imshow(CMapRot, "Parent", ax_2dmap);

end


function CMapShowButtonPushed(ax_2dmap)
%Displaying colormap in a separate window
CMapRot = getimage(ax_2dmap);
figure, imshow(CMapRot);
end


function SurfShowButtonPushed(nTerms, ax_2dmap)
global ns ms mags angles R Theta

%Compute cartesian coordinates
x = R.*cos(Theta);
y = R.*sin(Theta);

%Compute z values
z = zEval();

%Get colormap
colormap = getimage(ax_2dmap);
colormap = imrotate(colormap, 90);

%Find z limits
zMax = max(z, [], 'all');
zMin = -zMax;

%Plot surface
figure
surf(x,y,z,colormap, "EdgeColor", "none", "FaceColor", "interp");
zlim([zMin-5, zMax+5]);
axis manual off

%Show wave function
eq = "f(r,) = ";
for i = 1:nTerms
    newTerm = "";
    n = ns(i);
    m = ms(i);
    A = mags(i);
    phi = angles(i);
    newTerm = "J_" + n + "(_" + n + ",_"+ m+"r)"+"("+A+"cos(+"+phi+")"+"+sin(+"+phi+")"+")";
    eq = eq + newTerm;
    if i < nTerms
        eq = eq + "+";
    end
end
title(eq,"FontWeight","bold", "FontAngle",'italic','FontSize',9);

end


function animateWave(ax_2dmap)
global ns ms mags angles Num R Theta

%Calculate domain values
x = R.*cos(Theta);
y = R.*sin(Theta);

%Define wave function
z = @(t) zeros(Num+1, Num+1);
nTerms = numel(ns);
for k = 1:nTerms
    n = ns(k);
    m = ms(k);
    A = mags(k);
    phi = angles(k);
    L = besselzero(n,20,1);
    z = @(t) z(t) + besselj(n, L(m).*R).*(A.*(cos(n.*(Theta+phi)) + sin(n.*(Theta+phi)))).*cos(L(m).*t);
end

%Get colormap
colormap = getimage(ax_2dmap);

%Creates data structure for time series
step = 0.01;
upper = 10;
ts = 1:step:upper;
loops = numel(ts);

%Initialize video object
myVideo = VideoWriter("vibrating_drumhead"); 
myVideo.FrameRate = 15;
open(myVideo)

%Find limits for z axes.
zmax = max(zEval(), [], 'all');
zmin = -zmax;

figure
%Creates frames
for j = 1:loops
    %Plots surface on current frame and stores the frame to video object
    t = ts(j);
    surf(x, y, z(t), colormap, "EdgeColor","interp","FaceColor", "interp");
    zlim([zmin-5 zmax+5]);
    axis manual off
    frame = getframe(gcf);
    writeVideo(myVideo, frame);
end

close(myVideo);

end


function Y = scale(X, new_min, new_max)
%Helper function for scaling
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
