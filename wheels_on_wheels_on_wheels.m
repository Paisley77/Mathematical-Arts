%Constant variable
t = 0:0.001:2*pi;
%%
%User-defined Application Parameters
pFold = input("How many folds? ");
mType = input("Of what type? ");
nTerms = input("number of terms: ");
%%
%Set up window

%Set up figure
fig = uifigure("WindowStyle","normal","Name","CreateSymmetryUI","Position",[100 100 600 500]);

%Customize figure layout
g = uigridlayout(fig,[nTerms*2+1,3]);
g.RowHeight(1) = {20};
rowHeight = {480/(2*nTerms)};
g.RowHeight(2:nTerms*2+1) = repmat(rowHeight,1,nTerms*2);
g.ColumnWidth = {100,100,'3x'};

%Set up axis for the graph
ax = uiaxes(g);
ax.Layout.Row = [2 nTerms*2+1];
ax.Layout.Column = 3;
%%
%Set labels

%coef
coefLabel = uilabel(g,"Text","coef.(real&imag)");
coefLabel.HorizontalAlignment = 'center';
coefLabel.Layout.Row = 1;
coefLabel.Layout.Column = 1;
coefLabel.FontName = "TimesNewRoman";
coefLabel.FontSize = 10;

%freq
freqLabel = uilabel(g,"Text","freq.("+mType+" mod "+pFold+")");
freqLabel.HorizontalAlignment = 'center';
freqLabel.Layout.Row = 1;
freqLabel.Layout.Column = 2;
freqLabel.FontName = "TimesNewRoman";
freqLabel.FontSize = 10;

%graph
graphLabel = uilabel(g,"Text","Figure ("+pFold+" fold symmetry of type "+mType+")");
graphLabel.HorizontalAlignment = 'center';
graphLabel.Layout.Row = 1;
graphLabel.Layout.Column = 3;
grahphLabel.FontName = "TimesNewRoman";
graphLabel.FontSize = 10;
%%
%Initialize Figure
setFreq(nTerms,pFold,mType);
setCoef(nTerms);
createSliders(g,nTerms,pFold,mType,t,ax);
drawGraph(t,ax)
%% 
% Function Definitions

function setFreq(nTerms,pFold,mType)
    global freq
    freq = zeros(nTerms);
    for i = 1:nTerms
        freq(i) = mType + (i-1)*pFold;
    end
end

function setCoef(nTerms)
    global Rcoef Icoef coef
    Rcoef = ones(nTerms);
    Icoef = ones(nTerms);
    coef = complex(Rcoef, Icoef);
end


function createSliders(g,nTerms,pFold,mType,t,ax)
    for i = 1:nTerms
        %coef sliders
        %real part
        sld = uislider(g,"Value",1, "Limits", [-25 25], "ValueChangingFcn", ...
            @(src,event)updateRealCoef(src,event,i,t,ax));
        sld.Layout.Row = 2*i;
        sld.Layout.Column = 1;
        %imaginary part
        sld = uislider(g,"Value",1,"Limits", [-25 25], "ValueChangingFcn",...
            @(src,event)updateImagCoef(src,event,i,t,ax));
        sld.Layout.Row = 2*i+1;
        sld.Layout.Column = 1;
        %freq sliders
        sld = uislider(g,"Value",mType+(i-1)*pFold, "Limits", [-100 100], "ValueChangingFcn", ...
            @(src,event)updateFreq(src,event,pFold, mType, i,t,ax));
        sld.Layout.Row = [2*i 2*i+1];
        sld.Layout.Column = 2;
    end
end


function updateRealCoef(src,event,i,t,ax)
    global Rcoef Icoef coef
    Rcoef(i) = event.Value;
    coef = complex(Rcoef, Icoef);
    drawGraph(t,ax)
end


function updateImagCoef(src,event,i,t,ax)
    global Rcoef Icoef coef
    Icoef(i) = event.Value;
    coef = complex(Rcoef, Icoef);
    drawGraph(t,ax)
end


function updateFreq(src,event,pFold, mType, i,t,ax)
    global freq
    intVal = round(event.Value);
    %Set the freq. value to the nearest k = m mod p
    freq(i) = intVal - (mod(intVal,pFold) - mType);
    drawGraph(t,ax);
end


function drawGraph(t,ax)
    %Plot the graph
    global freq
    global coef
    X_coords = zeros(size(t));
    Y_coords = zeros(size(t));
    for i = 1:size(t,2)
        [X,Y] = sumOverExpr(t(i), freq, coef);
        X_coords(i) = X;
        Y_coords(i) = Y;
    end
    plot(ax,X_coords, Y_coords, "LineWidth", 1.5, "Color", [195/255, 50/255, 186/255])
    axis equal
    %display equation
    eq = "";
    for i = 1:size(freq,2)
        exponent = ""+freq(i)+"t";
        if i>1
            eq = eq + '+';
        end
        eq = eq + "(" + coef(i) + ")" + 'e' + exponent;
    end
    text(ax,ax.XLim(1),ax.YLim(2)-0.2,eq, 'HorizontalAlignment', 'left');
end


function [x,y] = extractXY(a,n)
    inner = complex(cos(n), sin(n));
    x = real(a * inner);
    y = imag(a * inner);
end


function [X,Y] = sumOverExpr(t, freq, coef)
    X = 0;
    Y = 0;
    for i = 1:size(freq,2)
        [x,y] = extractXY(coef(i), freq(i)*t);
        X = X + x;
        Y = Y + y;
    end
end