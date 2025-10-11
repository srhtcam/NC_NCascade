% NC^4_4 MATLAB surface patch Generator
% Disclaimer: This script generates only NC4 surface patches corresponding to n = 4 narrowing cascades.
% For other narrowing-cascade configurations, please refer to the Rhino–Grasshopper implementation.

clc;
clear all;

syms u v

% Read the CSV file (each row = [x y z])
[file, path] = uigetfile('*.txt', 'Select the .txt file which contain control points');
pointss = readmatrix(fullfile(path, file));
[~, name, ~] = fileparts(file);   % name = 'aaa'
outputFilepath = fullfile(path, name + ".bv");

% Check number of points
if size(pointss, 1) ~= 52
    error("The file must contain exactly 52 points. Current file contains %d points.", size(pointss,1));
end
disp('Control points loaded successfully:');

% Extended net surf creation
netpts = NetPtsOut(pointss);
tensorpoints = SetExtendedSurf(pointss);
Setbvfile(tensorpoints,outputFilepath);


% Number of points in file
numPts = size(netpts, 1);

nn = 6;
ii = 1;
mm = 0;
for i = 1:nn
    for j = 1:(nn - mm) 
        if ii <= size(netpts,1)
            d_ij{j, i} = netpts(ii, :);  % store [x y z] vector
            ii = ii + 1;
        else
            d_ij{j, i} = [NaN NaN NaN];
        end
    end
    for j = (nn - mm + 1):nn % Fill the rest with empty (unset) entries
        d_ij{j, i} = [NaN NaN NaN];
    end
    if i > 1 && i < nn - 1
        mm = mm + 1;
    end
end


% NC4 (n = 4) cap surf creation
tij = cell(9, 9);
[t_ij, tT_ij, tB_ij] = SetTensorBorder(d_ij);
[rp_t_ij, rp_tB_ij, rp_tT_ij] = SetReparametizedTensor(t_ij, tT_ij, tB_ij);
pts = RepTensorPtOut(rp_t_ij, rp_tB_ij, rp_tT_ij);


m = 15; n = 9;
CP = FillControlPoints(pts);



% % % % % % % 
% Calculation of c_ij points
% 
% Define the matrix
M = (1/128) * [ ...
    4  88 36  0  0  0  0;
    0  36 88  4  0  0  0;
    0   0  0  0 49 78  1;
    0   0  0  0  1 78  49];

% Define the vector of d-values
dd = [d_ij(1,4); d_ij(2,4); d_ij(3,4); d_ij(4,4); d_ij(2,3); d_ij(3,3); d_ij(4,3)];

% Initialize output c vectors as 1x3 zero vectors
c = cell(4,1);
for i = 1:4
    c{i} = zeros(1,3);  % [x y z]
    for j = 1:7
        c{i} = c{i} + M(i,j) * dd{j}; 
    end
end

% Assign to individual variables
c23 = cell2mat(c(1));
c33 = cell2mat(c(2));
c22 = cell2mat(c(3));
c32 = cell2mat(c(4));
% % % % % % % 


% Defining contiunity rules
CP = ContiunityRules(CP,c23,c33,c22,c32);

% all point list
ptlist = [];
for ifor=1:m
    for jfor = 1:n
        ptlist = [ptlist; CP{ifor,jfor}];
    end
end
% 


NC4Surface(ptlist,outputFilepath);

% end of the code








% % % % % % % % % % 
% Functions       % 
% % % % % % % % % % 


% % % % % % % % % % % % % % % % % % % % % % % % % 
% Contiunity rules....
% 
function Rules = ContiunityRules(CP,MTT4,MTT6,MT4,MT6)

CP{8,4} = 0.5.*(MTT4 + MT4);
CP{8,6} = 0.5.*(MTT6 + MT6);

% diamond rule
% top part bi-2
A1_4 = (2.*CP{14,4}) - CP{15,4};
CP{11,4} = 0.5.*(A1_4 + MTT4);
CP{12,4} = 0.5.*(A1_4 + CP{11,4});
CP{13,4} = (1/6).*CP{15,4} + (2/3).*A1_4 + (1/6).*CP{11,4};

A1_6 = (2.*CP{14,6}) - CP{15,6};
CP{11,6} = 0.5.*(A1_6 + MTT6);
CP{12,6} = 0.5.*(A1_6 + CP{11,6});
CP{13,6} = (1/6).*CP{15,6} + (2/3).*A1_6 + (1/6).*CP{11,6};

% bottom part bi-2
B1_4 = (2.*CP{2,4}) - CP{1,4};
CP{5,4} = 0.5.*(B1_4 + MT4);
CP{4,4} = 0.5.*(B1_4 + CP{5,4});
CP{3,4} = (1/6).*CP{1,4} + (2/3).*B1_4 + (1/6).*CP{5,4};

B1_6 = (2.*CP{2,6}) - CP{1,6};
CP{5,6} = 0.5.*(B1_6 + MT6);
CP{4,6} = 0.5.*(B1_6 + CP{5,6});
CP{3,6} = (1/6).*CP{1,6} + (2/3).*B1_6 + (1/6).*CP{5,6};

% diamond rule
% 

% top spine points
CP{9,4}  = (1/3).*CP{8,4} + (2/3).*MTT4;
CP{10,4}  = (2/3).*MTT4 + (1/3).*CP{11,4};

CP{9,6}  = (1/3).*CP{8,6} + (2/3).*MTT6;
CP{10,6}  = (2/3).*MTT6 + (1/3).*CP{11,6};

%mid spine points
CP{6,4}  = (1/3).*CP{5,4} + (2/3).*MT4;
CP{7,4}  = (2/3).*MT4 + (1/3).*CP{8,4};

CP{6,6}  = (1/3).*CP{5,6} + (2/3).*MT6;
CP{7,6}  = (2/3).*MT6 + (1/3).*CP{8,6};

% g1 equations horzontal
CP{3,3}  =  0.5.*(CP{3,2} + CP{3,4});
CP{4,3}  =  0.5.*(CP{4,2} + CP{4,4});
CP{5,3}  =  0.5.*(CP{5,2} + CP{5,4});
CP{6,3}  =  0.5.*(CP{6,2} + CP{6,4});
CP{7,3}  =  0.5.*(CP{7,2} + CP{7,4});
CP{8,3}  =  0.5.*(CP{8,2} + CP{8,4});
CP{9,3}  =  0.5.*(CP{9,2} + CP{9,4});
CP{10,3} =  0.5.*(CP{10,2} + CP{10,4});
CP{11,3} =  0.5.*(CP{11,2} + CP{11,4});
CP{12,3} =  0.5.*(CP{12,2} + CP{12,4});
CP{13,3} =  0.5.*(CP{13,2} + CP{13,4});

CP{3,5}  =  0.5.*(CP{3,4} + CP{3,6});
CP{4,5}  =  0.5.*(CP{4,4} + CP{4,6});
CP{5,5}  =  0.5.*(CP{5,4} + CP{5,6});
CP{6,5}  =  0.5.*(CP{6,4} + CP{6,6});
CP{7,5}  =  0.5.*(CP{7,4} + CP{7,6});
CP{8,5}  =  0.5.*(CP{8,4} + CP{8,6});
CP{9,5}  =  0.5.*(CP{9,4} + CP{9,6});
CP{10,5} =  0.5.*(CP{10,4} + CP{10,6});
CP{11,5} =  0.5.*(CP{11,4} + CP{11,6});
CP{12,5} =  0.5.*(CP{12,4} + CP{12,6});
CP{13,5} =  0.5.*(CP{13,4} + CP{13,6});

CP{3,7}  =  0.5.*(CP{3,6} + CP{3,8});
CP{4,7}  =  0.5.*(CP{4,6} + CP{4,8});
CP{5,7}  =  0.5.*(CP{5,6} + CP{5,8});
CP{6,7}  =  0.5.*(CP{6,6} + CP{6,8});
CP{7,7}  =  0.5.*(CP{7,6} + CP{7,8});
CP{8,7}  =  0.5.*(CP{8,6} + CP{8,8});
CP{9,7}  =  0.5.*(CP{9,6} + CP{9,8});
CP{10,7} =  0.5.*(CP{10,6} + CP{10,8});
CP{11,7} =  0.5.*(CP{11,6} + CP{11,8});
CP{12,7} =  0.5.*(CP{12,6} + CP{12,8});
CP{13,7} =  0.5.*(CP{13,6} + CP{13,8});
% 

Rules = CP;
end
% % % % % % % % % % % % % % % % % % % % % % % 


% % % 
% Defined functions for this matlab code
% 

function [newControlPtUpper, newControlPtBottom] = DecasteljauDivision(controlpt, u)
    % controlpt: Nx3 array of control points [x y z]
    % u: parameter in [0,1]
    % newControlPtUpper, newControlPtBottom: Nx3 arrays (the two halves)

    n = size(controlpt, 1);
    
    % Initialize output arrays
    newControlPtUpper = cell(n, 1);
    newControlPtBottom = cell(n, 1);
    
    % Initialize point matrix (cell array of levels)
    ptMatrix = cell(n, 1);
    ptMatrix{1} = controlpt;
    
    % De Casteljau iteration
    for i = 1:n
        % Create next level
        ptMatrix{i+1} = cell(size(ptMatrix{i}, 1) - 1, 1);
        for j = 1:(size(ptMatrix{i}, 1) - 1)
            ptMatrix{i+1}{j, :} = (1 - u) * ptMatrix{i}{j, :} + u * ptMatrix{i}{j + 1, :};
        end
        
        % Assign upper and bottom control points
        newControlPtUpper{i, :} = ptMatrix{i}{1, :};
        newControlPtBottom{n - i + 1, :} = ptMatrix{i}{end, :};
    end
end


function [t_ij, tT_ij, tB_ij] = SetTensorBorder(d_ij)
% SetTensorBorder constructs tensor border structures from control net d_ij
% d_ij: {m x n} cell array, each cell is [x y z]
% Output:
%   t_ij  - {1xK} cell array of 2x3x3 numeric blocks
%   tT_ij - {1x1} top tensor block
%   tB_ij - {1x4} bottom tensor blocks

t_ij = {};
tT_ij = {};
tB_ij = {};

% Helper inline functions
SetQMiddle = @(p0,p1,p2,p3) (p0 + p1 + p2 + p3)/4.0;
SetEMiddle = @(p0,p1) (p0 + p1)/2.0;

% LEFT BORDER
for i = 1:4
    t = cell(2, 3);
    t{1,1} = SetQMiddle(d_ij{1,i}, d_ij{1,i+1}, d_ij{2,i}, d_ij{2,i+1});
    t{2,1} = SetEMiddle(d_ij{2,i}, d_ij{2,i+1});
    t{1,2} = SetEMiddle(d_ij{1,i+1}, d_ij{2,i+1});
    t{2,2} = d_ij{2,i+1};
    t{1,3} = SetQMiddle(d_ij{1,i+1}, d_ij{1,i+2}, d_ij{2,i+1}, d_ij{2,i+2});
    t{2,3} = SetEMiddle(d_ij{2,i+1}, d_ij{2,i+2});
    t_ij{end+1} = t;
end

% RIGHT BORDER
t = cell(2,3);
t{1,1} = SetEMiddle(d_ij{5,1}, d_ij{5,2});
t{2,1} = SetQMiddle(d_ij{5,1}, d_ij{5,2}, d_ij{6,1}, d_ij{6,2});
t{1,2} = d_ij{5,2};
t{2,2} = SetEMiddle(d_ij{5,2}, d_ij{6,2});
t{1,3} = SetEMiddle(d_ij{5,2}, d_ij{4,3});
t{2,3} = SetQMiddle(d_ij{5,2}, d_ij{4,3}, d_ij{6,2}, d_ij{5,3});
t_ij{end+1} = t;

t = cell(2,3);
t{1,1} = SetEMiddle(d_ij{5,2}, d_ij{4,3});
t{2,1} = SetQMiddle(d_ij{5,2}, d_ij{4,3}, d_ij{6,2}, d_ij{5,3});
t{1,2} = d_ij{4,3};
t{2,2} = SetEMiddle(d_ij{4,3}, d_ij{5,3});
t{1,3} = SetEMiddle(d_ij{4,3}, d_ij{3,4});
t{2,3} = SetQMiddle(d_ij{4,3}, d_ij{3,4}, d_ij{5,3}, d_ij{4,4});
t_ij{end+1} = t;

t = cell(2,3);
t{1,1} = SetEMiddle(d_ij{4,3}, d_ij{3,4});
t{2,1} = SetQMiddle(d_ij{4,3}, d_ij{3,4}, d_ij{5,3}, d_ij{4,4});
t{1,2} = d_ij{3,4};
t{2,2} = SetEMiddle(d_ij{3,4}, d_ij{4,4});
t{1,3} = SetEMiddle(d_ij{3,4}, d_ij{2,5});
t{2,3} = SetQMiddle(d_ij{3,4}, d_ij{2,5}, d_ij{4,4}, d_ij{3,5});
t_ij{end+1} = t;

t = cell(2,3);
t{1,1} = SetEMiddle(d_ij{3,4}, d_ij{2,5});
t{2,1} = SetQMiddle(d_ij{3,4}, d_ij{2,5}, d_ij{4,4}, d_ij{3,5});
t{1,2} = d_ij{2,5};
t{2,2} = SetEMiddle(d_ij{2,5}, d_ij{3,5});
t{1,3} = SetEMiddle(d_ij{2,5}, d_ij{2,6});
t{2,3} = SetQMiddle(d_ij{2,5}, d_ij{2,6}, d_ij{3,5}, d_ij{3,6});
t_ij{end+1} = t;

% TOP BORDER
t = cell(3,3);
t{1,1} = SetEMiddle(d_ij{1,5}, d_ij{2,5});
t{2,1} = d_ij{2,5};
t{3,1} = SetEMiddle(d_ij{2,5}, d_ij{3,5});
t{1,2} = SetQMiddle(d_ij{1,5}, d_ij{1,6}, d_ij{2,5}, d_ij{2,6});
t{2,2} = SetEMiddle(d_ij{2,5}, d_ij{2,6});
t{3,2} = SetQMiddle(d_ij{2,5}, d_ij{2,6}, d_ij{3,5}, d_ij{3,6});
tT_ij{1} = t;

% BOTTOM BORDER
for i = 1:4
    t = cell(3,2);
    t{1,1} = SetQMiddle(d_ij{i,1}, d_ij{i,2}, d_ij{i+1,1}, d_ij{i+1,2});
    t{2,1} = SetEMiddle(d_ij{i+1,1}, d_ij{i+1,2});
    t{3,1} = SetQMiddle(d_ij{i+1,1}, d_ij{i+1,2}, d_ij{i+2,1}, d_ij{i+2,2});
    t{1,2} = SetEMiddle(d_ij{i,2}, d_ij{i+1,2});
    t{2,2} = d_ij{i+1,2};
    t{3,2} = SetEMiddle(d_ij{i+1,2}, d_ij{i+2,2});
    tB_ij{end+1} = t;
end

end


function [rp_t_ij, rp_tB_ij, rp_tT_ij] = SetReparametizedTensor(t_ij, tT_ij, tB_ij)
% Reparametrize tensor borders for left/right/top/bottom patches
% Input:
%   t_ij, tT_ij, tB_ij - cell arrays containing control point blocks (each 2x3 cell of [x y z])
% Output:
%   rp_t_ij, rp_tB_ij, rp_tT_ij - reparameterized cell arrays

% Initialize outputs
rp_t_ij = cell(1, 8);
rp_tB_ij = cell(1, 4);
rp_tT_ij = cell(1, 4);


% Define reparameterization coefficients
a0 = [1.0];
a1 = [1.0];
a2 = [1.0 - (1.0 / (2*4))];
for i = 2:4-1
    a0 = [a0, 1.0 - (((2*i) - 3.0) / (2*4))];
    a1 = [a1, 1.0 - (((2*(i+1)) - 3.0) / (2*4))];
    a2 = [a2, 0];
end
a0 = [a0, 1.0 - (((2*4) - 3.0) / (2*4))];
a1 = [a1, 1.0 / 4];
a2 = [a2, 1.0 / 4];

% Duplicate for right tensors (8 total)
a0 = [a0, a0];
a1 = [a1, a1];
a2 = [a2, a2];

% LEFT + RIGHT TENSORS
for i = 1:8
    % Extract numeric [x y z] points
    t00 = t_ij{i}{1,1}; 
    t10 = t_ij{i}{1,2}; 
    t20 = t_ij{i}{1,3};
    t01 = t_ij{i}{2,1}; 
    t11 = t_ij{i}{2,2}; 
    t21 = t_ij{i}{2,3};

    % For right tensors
    if i > 4
        t00 = t_ij{i}{2,1}; 
        t10 = t_ij{i}{2,2}; 
        t20 = t_ij{i}{2,3};
        t01 = t_ij{i}{1,1}; 
        t11 = t_ij{i}{1,2}; 
        t21 = t_ij{i}{1,3};
    end

    if ismember(i, [1, 4, 5, 8])
        par_t00 = t00;
        par_t10 = (t00 + t10) / 2;
        par_t20 = (t00 / 6) + (2 * t10 / 3) + (t20 / 6);
        par_t30 = (t10 + t20) / 2;
        par_t40 = t20;

        par_t01 = -a0(i)*t00 + a0(i)*t01 + t00;
        par_t11 = -a0(i)*t10/2 + a0(i)*t11/2 - a1(i)*t00/2 + a1(i)*t01/2 + (t00 + t10)/2;
        par_t21 = -a0(i)*t20/6 + a0(i)*t21/6 - 2*a1(i)*t10/3 + 2*a1(i)*t11/3 - a2(i)*t00/6 + a2(i)*t01/6 + (t00/6 + 2*t10/3 + t20/6);
        par_t31 = -a1(i)*t20/2 + a1(i)*t21/2 - a2(i)*t10/2 + a2(i)*t11/2 + (t10 + t20)/2;
        par_t41 = -a2(i)*t20 + a2(i)*t21 + t20;

        rp_t_ij{i} = {par_t00, par_t10, par_t20, par_t30, par_t40; ...
                      par_t01, par_t11, par_t21, par_t31, par_t41};
    else
        par_t00 = t00;
        par_t10 = (t00/3) + (2*t10/3);
        par_t20 = (2*t10/3) + (t20/3);
        par_t30 = t20;

        par_t01 = -a0(i)*t00 + a0(i)*t01 + t00;
        par_t11 = -2*a0(i)*t10/3 + 2*a0(i)*t11/3 - a1(i)*t00/3 + a1(i)*t01/3 + (t00/3 + 2*t10/3);
        par_t21 = -a0(i)*t20/3 + a0(i)*t21/3 - 2*a1(i)*t10/3 + 2*a1(i)*t11/3 + (2*t10/3 + t20/3);
        par_t31 = -a1(i)*t20 + a1(i)*t21 + t20;

        rp_t_ij{i} = {par_t00, par_t10, par_t20, par_t30; ...
                      par_t01, par_t11, par_t21, par_t31};
    end
end

% BOTTOM TENSORS
for i = 1:4
    tB = tB_ij{i};
    rp_tB_ij{i} = cell(3,2);
    for r = 1:3
        rp_tB_ij{i}{r,1} = tB{r,1};
        rp_tB_ij{i}{r,2} = (tB{r,1} + tB{r,2}) / 2;
    end
end

%  TOP TENSOR
% bottom and top decasteljau subdivision
newcntrolpts_bottom = {};
newcntrolpts_top = {};
decastel_bottom = { [tT_ij{1}(1,1,:)' ; tT_ij{1}(2,1,:)' ; tT_ij{1}(3,1,:)'] };
decastel_top    = { [tT_ij{1}(1,2,:)' ; tT_ij{1}(2,2,:)' ; tT_ij{1}(3,2,:)'] };

for i = 4:-1:2
    [new0_b, new1_b] = DecasteljauDivision(decastel_bottom{end}, 1.0/i);
    decastel_bottom{end+1} = new1_b;
    newcntrolpts_bottom{end+1} = new0_b;

    [new0_t, new1_t] = DecasteljauDivision(decastel_top{end}, 1.0/i);
    decastel_top{end+1} = new1_t;
    newcntrolpts_top{end+1} = new0_t;
end

newcntrolpts_bottom{end+1} = decastel_bottom{end};
newcntrolpts_top{end+1} = decastel_top{end};

for i = 1:4
    rp_tT_ij{i} = cell(3,2);
    for k = 1:3
        % Extract numeric rows from bottom and top
        bottom_val = newcntrolpts_bottom{i}{k};  % assume each cell contains numeric row [x y z]
        top_val    = newcntrolpts_top{i}{k};

        % Compute average
        rp_tT_ij{i}{k,1} = (bottom_val + top_val)/2;
        rp_tT_ij{i}{k,2} = top_val;
    end
end
end


function pts = RepTensorPtOut(rp_t_ij, rp_tB_ij, rp_tT_ij)
% Inputs:
%   rp_t_ij, rp_tB_ij, rp_tT_ij : cell arrays, each cell contains a [m×n×3] numeric array
% Output:
%   pts : [N×3] numeric matrix of concatenated control points

pts = [];

% vertical line 0
for i = 1:4
    if i == 1
        pts = [pts; squeeze(rp_t_ij{i}(1,1,:))'];
    end
    pts = [pts; squeeze(rp_t_ij{i}(1,2,:))'];
    pts = [pts; squeeze(rp_t_ij{i}(1,3,:))'];
    pts = [pts; squeeze(rp_t_ij{i}(1,4,:))'];
    if i == 1 || i == 4
        pts = [pts; squeeze(rp_t_ij{i}(1,5,:))'];
    end
end

% vertical line 1
for i = 1:4
    if i == 1
        pts = [pts; squeeze(rp_t_ij{i}(2,1,:))'];
    end
    pts = [pts; squeeze(rp_t_ij{i}(2,2,:))'];
    pts = [pts; squeeze(rp_t_ij{i}(2,3,:))'];
    pts = [pts; squeeze(rp_t_ij{i}(2,4,:))'];
    if i == 1 || i == 4
        pts = [pts; squeeze(rp_t_ij{i}(2,5,:))'];
    end
end

% mixed B and T tensors
for i = 1:3
    if i ~= 1
        pts = [pts; squeeze(rp_tB_ij{i}(2,1,:))'];
        pts = [pts; squeeze(rp_tB_ij{i}(2,2,:))'];
        pts = [pts; squeeze(rp_tT_ij{i}(2,1,:))'];
        pts = [pts; squeeze(rp_tT_ij{i}(2,2,:))'];
    end
    if i ~= 4
        pts = [pts; squeeze(rp_tB_ij{i}(3,1,:))'];
        pts = [pts; squeeze(rp_tB_ij{i}(3,2,:))'];
        pts = [pts; squeeze(rp_tT_ij{i}(3,1,:))'];
        pts = [pts; squeeze(rp_tT_ij{i}(3,2,:))'];
    end
end

% vertical line 2n-1
for i = 5:8
    if i == 5
        pts = [pts; squeeze(rp_t_ij{i}(2,1,:))'];
    end
    pts = [pts; squeeze(rp_t_ij{i}(2,2,:))'];
    pts = [pts; squeeze(rp_t_ij{i}(2,3,:))'];
    pts = [pts; squeeze(rp_t_ij{i}(2,4,:))'];
    if i == 5 || i == 8
        pts = [pts; squeeze(rp_t_ij{i}(2,5,:))'];
    end
end

% vertical line 2n
for i = 5:8
    if i == 5
        pts = [pts; squeeze(rp_t_ij{i}(1,1,:))'];
    end
    pts = [pts; squeeze(rp_t_ij{i}(1,2,:))'];
    pts = [pts; squeeze(rp_t_ij{i}(1,3,:))'];
    pts = [pts; squeeze(rp_t_ij{i}(1,4,:))'];
    if i == 5 || i == 8
        pts = [pts; squeeze(rp_t_ij{i}(1,5,:))'];
    end
end

end


function CP = FillControlPoints(pts)
% pts : [N×3] numeric matrix (like output of RepTensorPtOut)
% CP  : {15×9} cell array, where each CP{i,j} = [x y z]

% Create a symbolic cell array structure CP
CP = cell(15, 9);

% Define symbolic variables for each element in CP
for i = 1:15
    for j = 1:9
        CP{i, j} = sym(['P' num2str(i) num2str(j)]);
    end
end

% helper to convert from 0-based (C#) to 1-based (MATLAB)
toMat = @(idx) idx + 1;

% fill first column (1)
for k = 0:14
    CP{k+1,1} = cell2mat(pts(toMat(k), :));
end

% fill second column (2)
for k = 15:29
    CP{k-14,2} = cell2mat(pts(toMat(k), :));
end

% column 3
CP{1 ,3}  = cell2mat(pts(toMat(30), :));
CP{2 ,3}  = cell2mat(pts(toMat(31), :));
CP{14,3}  = cell2mat(pts(toMat(32), :));
CP{15,3}  = cell2mat(pts(toMat(33), :));

% column 4
CP{1 ,4}  = cell2mat(pts(toMat(34), :));
CP{2 ,4}  = cell2mat(pts(toMat(35), :));
CP{14,4}  = cell2mat(pts(toMat(36), :));
CP{15,4}  = cell2mat(pts(toMat(37), :));

% column 5
CP{1 ,5}  = cell2mat(pts(toMat(38), :));
CP{2 ,5}  = cell2mat(pts(toMat(39), :));
CP{14,5}  = cell2mat(pts(toMat(40), :));
CP{15,5}  = cell2mat(pts(toMat(41), :));

% column 6
CP{1 ,6}  = cell2mat(pts(toMat(42), :));
CP{2 ,6}  = cell2mat(pts(toMat(43), :));
CP{14,6}  = cell2mat(pts(toMat(44), :));
CP{15,6}  = cell2mat(pts(toMat(45), :));

% column 7
CP{1 ,7}  = cell2mat(pts(toMat(46), :));
CP{2 ,7}  = cell2mat(pts(toMat(47), :));
CP{14,7}  = cell2mat(pts(toMat(48), :));
CP{15,7}  = cell2mat(pts(toMat(49), :));

% column 8
for k = 50:64
    CP{k-49,8} = cell2mat(pts(toMat(k), :));
end

% column 9
for k = 65:79
    CP{k-64,9} = cell2mat(pts(toMat(k), :));
end
end


function NC4Surface(points,path)
    % points: Nx3 matrix [X Y Z]
    % Output: writes text to NC4.bv and returns bvfile string

    bvfile = "Group 2 NC4cap\n";

    %Helper function for formatted points
    function str = appendPoints(idx, pts)
        str = join(string(pts(idx,1)) + " " + string(pts(idx,2)) + " " + string(pts(idx,3)) + "\n", "");
    end

    % Loop 1: ii = 0:2:6
    for ii = 1:2:7
        i = ii;
        bvfile = bvfile + "5\n";
        bvfile = bvfile + "4 2\n";
        idx = [i i+1 i+2 i+9 i+10 i+11 i+18 i+19 i+20 ...
                i+27 i+28 i+29 i+36 i+37 i+38];
        bvfile = bvfile + appendPoints(idx,points);
    end

    % Loop 2: ii = 36:2:42
    for ii = 37:2:43
        i = ii;
        bvfile = bvfile + "5\n";
        bvfile = bvfile + "3 2\n";
        idx = [i i+1 i+2 i+9 i+10 i+11 i+18 i+19 i+20 ...
                i+27 i+28 i+29];
        bvfile = bvfile + appendPoints(idx,points);
    end

    % Loop 3: ii = 63:2:69
    for ii = 64:2:70
        i = ii;
        bvfile = bvfile + "5\n";
        bvfile = bvfile + "3 2\n";
        idx = [i i+1 i+2 i+9 i+10 i+11 i+18 i+19 i+20 ...
                i+27 i+28 i+29];
        bvfile = bvfile + appendPoints(idx,points);
    end

    % Loop 4: ii = 90:2:96
    for ii = 91:2:97
        i = ii;
        bvfile = bvfile + "5\n";
        bvfile = bvfile + "4 2\n";
        idx = [i i+1 i+2 i+9 i+10 i+11 i+18 i+19 i+20 ...
                i+27 i+28 i+29 i+36 i+37 i+38];
        bvfile = bvfile + appendPoints(idx,points);
    end

    %Write to .bv file
    fid = fopen(path, 'a');
    if fid == -1
        error('Cannot open file for appending: %s', path +"FC4SurfaceOutput.bv");
    end
    fprintf(fid, bvfile);
    fclose(fid);
end


function netpts = NetPtsOut(net_extent)
    netpts = [];
    
    i = 10;
    ii = 0;
    iii = 0;
    
    while i < (size(net_extent, 1) - 5)

        for j = 1:(6 - iii)
            netpts = [netpts; net_extent(i, :)];
            i = (i + 1);
        end
        
        if ii > 0 && ii < (4)
            iii = iii + 1;
        end
        ii = ii + 1;
        i = (i + 2);
    end
end


function p = SetEMiddle(p1, p2)
    p = (p1 + p2) / 2;
end


function p = SetQMiddle(p1, p2, p3, p4)
    p = (p1 + p2 + p3 + p4) / 4;
end


function tensorpoints = SetExtendedSurf(net_extent)

    tensorpoints = [];

    for ii = 0:1
        if ii ~= 0
            for i = (ii * 8):(ii * 8) + 6
                tensorpoints = [tensorpoints; SetEMiddle(net_extent(i + 1, :), net_extent(i + 2, :))];
                if i ~= (ii * 8) + 6
                    tensorpoints = [tensorpoints; net_extent(i + 2, :)];
                end
            end
        end

        if ii == 7
            continue;
        end

        for i = (ii * 8):(ii * 8) + 6
            if i ~= (ii * 8)
                tensorpoints = [tensorpoints; SetEMiddle(net_extent(i + 1, :), net_extent(i + 9, :))];
            end
            tensorpoints = [tensorpoints; SetQMiddle(net_extent(i + 1, :), net_extent(i + 2, :), net_extent(i + 9, :), net_extent(i + 10, :))];
        end
    end

    idx = [16];
    for ii = 1:length(idx)
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 1, :), net_extent(idx(ii) + 2, :))];
        tensorpoints = [tensorpoints; net_extent(idx(ii) + 2, :)];
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 3, :))];

        tensorpoints = [tensorpoints; SetQMiddle(net_extent(idx(ii) + 1, :), net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 9, :), net_extent(idx(ii) + 10, :))];
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 10, :))];
        tensorpoints = [tensorpoints; SetQMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 3, :), net_extent(idx(ii) + 10, :), net_extent(idx(ii) + 11, :))];
    end

    idx = [21, 24];
    for ii = 1:length(idx)
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 1, :), net_extent(idx(ii) + 2, :))];
        tensorpoints = [tensorpoints; net_extent(idx(ii) + 2, :)];
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 3, :))];

        tensorpoints = [tensorpoints; SetQMiddle(net_extent(idx(ii) + 1, :), net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 8, :), net_extent(idx(ii) + 9, :))];
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 9, :))];
        tensorpoints = [tensorpoints; SetQMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 3, :), net_extent(idx(ii) + 9, :), net_extent(idx(ii) + 10, :))];
    end

    idx = [28, 31];
    for ii = 1:length(idx)
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 1, :), net_extent(idx(ii) + 2, :))];
        tensorpoints = [tensorpoints; net_extent(idx(ii) + 2, :)];
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 3, :))];

        tensorpoints = [tensorpoints; SetQMiddle(net_extent(idx(ii) + 1, :), net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 7, :), net_extent(idx(ii) + 8, :))];
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 8, :))];
        tensorpoints = [tensorpoints; SetQMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 3, :), net_extent(idx(ii) + 8, :), net_extent(idx(ii) + 9, :))];
    end

    idx = [34, 37, 39];
    for ii = 1:length(idx)
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 1, :), net_extent(idx(ii) + 2, :))];
        tensorpoints = [tensorpoints; net_extent(idx(ii) + 2, :)];
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 3, :))];

        tensorpoints = [tensorpoints; SetQMiddle(net_extent(idx(ii) + 1, :), net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 6, :), net_extent(idx(ii) + 7, :))];
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 7, :))];
        tensorpoints = [tensorpoints; SetQMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 3, :), net_extent(idx(ii) + 7, :), net_extent(idx(ii) + 8, :))];

        if idx(ii) == 37
            tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 3, :), net_extent(idx(ii) + 8, :))];
        end
    end

    idx = [42, 44];
    for ii = 1:length(idx)
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 1, :), net_extent(idx(ii) + 2, :))];
        tensorpoints = [tensorpoints; net_extent(idx(ii) + 2, :)];
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 3, :))];

        if idx(ii) == 42
            tensorpoints = [tensorpoints; net_extent(idx(ii) + 3, :)];
        end

        tensorpoints = [tensorpoints; SetQMiddle(net_extent(idx(ii) + 1, :), net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 6, :), net_extent(idx(ii) + 7, :))];
        tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 7, :))];
        tensorpoints = [tensorpoints; SetQMiddle(net_extent(idx(ii) + 2, :), net_extent(idx(ii) + 3, :), net_extent(idx(ii) + 7, :), net_extent(idx(ii) + 8, :))];

        if idx(ii) == 42
            tensorpoints = [tensorpoints; SetEMiddle(net_extent(idx(ii) + 3, :), net_extent(idx(ii) + 8, :))];
        end
    end
end


function Setbvfile(tensorpoints,path)

    % tensorpoints: Nx3 matrix or cell array of [x y z] points
    if iscell(tensorpoints)
        tensorpoints = cell2mat(tensorpoints(:));
    end

    bvfile = {};
    bvfile{end+1} = 'Group 1 extend';
    nn = (3 + 3);
    nn2 = 2 * nn + 1;

    % Main surface tensor writing
    for ii = 0:(nn - 1)
        i = ii * 2;

        bvfile{end+1} = '5';
        bvfile{end+1} = '2 2';

        for k = [i, i+1, i+2, i+nn2, i+nn2+1, i+nn2+2, i+2*nn2, i+2*nn2+1, i+2*nn2+2]
            pt = tensorpoints(k+1, :);
            bvfile{end+1} = sprintf('%.6f %.6f %.6f', pt(1), pt(2), pt(3));
        end
    end

    % Build the index matrix
    index = cell(11, 1);
    index{1} = [26, 27, 28, 39, 40, 41, 42, 43, 44];
    index{2} = [36, 37, 38, 45, 46, 47, 48, 49, 50];
    index{3} = [42, 43, 44, 51, 52, 53, 54, 55, 56];
    index{4} = [48, 49, 50, 57, 58, 59, 60, 61, 62];
    index{5} = [54, 55, 56, 63, 64, 65, 66, 67, 68];
    index{6} = [60, 61, 62, 69, 70, 71, 72, 73, 74];
    index{7} = [66, 67, 68, 75, 76, 77, 78, 79, 80];
    index{8} = [72, 73, 74, 82, 83, 84, 85, 86, 87];
    index{9} = [78, 79, 80, 88, 89, 90, 92, 93, 94];
    index{10} = [80, 81, 85, 90, 91, 96, 94, 95, 99];
    index{11} = [85, 86, 87, 96, 97, 98, 99, 100, 101];


    % Write remaining tensor surfaces
    for i = 1:length(index)
        bvfile{end+1} = '5';
        bvfile{end+1} = '2 2';

        for k = 1:9
            idx = index{i}(k) + 1;
            pt = tensorpoints(idx, :);
            bvfile{end+1} = sprintf('%.6f %.6f %.6f', pt(1), pt(2), pt(3));
        end
    end

    % Write to .bv file
    fid = fopen(path, 'w');
    if fid == -1
        error('Failed to open or create file: %s', path +"FC4SurfaceOutput.bv");
    end
    for i = 1:length(bvfile)
        fprintf(fid, '%s\n', bvfile{i});
    end
    fclose(fid);
end
