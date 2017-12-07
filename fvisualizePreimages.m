% Discription:
% 
% fvisualizePreimages --> fuction that is used to visualize preimages.
%   Takes a variable number of arguments in the form of ('flag',input).
%   The number of preimages plotted is max([numel(theta),numel(phi)]).
%     If theta(preim) or phi(preim) exceeds the index, the last theta/phi
%     valus is used eg. theta(end) or phi(end).
%   Lighing is applied to each surface during creation.
%   Use ('clf',false) to use the figure without clearing (defalt is true).
% 
%     input: () if nargin = 0, nn from the base workspace is used
%     optional input: (nn, theta, phi, alpha, epsilon, fignum, clf)
%       nn: 4D double array with size [m,n,p,3]
%       theta: list of polar angles to use [-pi:pi]
%       phi: list of azamuthal angle to use [0:2*pi]
%       alpha: transparency of preimages [0:1]
%       epsilon: value used to generate isosurface in 'diffmag' [0:2]
%       fignum: figure to use [scalar integer from 1 to 2147483646]
%       clf: toggle clf true/false ([true],false,1,0)
%       az,el: view(az,el), default [0,90]
% 
% example of usage:
%       >> clear, nn = fknotansatz; fvisualizePreimages('alpha',0.7)
% 

function fvisualizePreimages(varargin)

% Defalt values:
nn = [];
theta = pi/2;
phi = linspace(0,2*pi-2*pi/5,5);
alpha = 1;
epsilon = 0.3;
fignum = 1;
clfIO = 1;
az = 0;
el = 90;

% Argument values
nVarargs = length(varargin);
for k = 1:2:nVarargs
    switch lower(varargin{k}) % case insensitive
        case 'nn'
            nn = varargin{k+1};
        case 'theta'
            theta = varargin{k+1};
        case 'phi'
            phi = varargin{k+1};
        case 'alpha'
            alpha = varargin{k+1};
        case 'epsilon'
            epsilon = varargin{k+1};
        case 'fignum'
            fignum = varargin{k+1};
        case 'clf'
            clfIO = varargin{k+1};
        case 'az'
            az = varargin{k+1};
        case 'el'
            el = varargin{k+1};
    end
end

% if nn is not given, look for nn in the base workspace
if nargin==0 || isempty(nn)
    try
        disp('Looking for ''nn'' in the base Workspace')
        nn = evalin('base','nn');
    catch
        error('no ''nn'' variable found in the base Workspace!')
    end
end

% if nn is wrong shape
[m,n,p,test] = size(nn);
if ~(test==3)
    error('First argument must be a 4D double array with size [m,n,p,3]!')
end

% ignore out of range theta/phi values, return if no valid value
k = 0;
while k<numel(theta)
    k = k+1;
    if theta(k) < -pi || theta(k) > pi
        warning('Out of range theta value ignored!')
        theta = theta(theta~=theta(k));
        if isempty(theta),error('No valid theta value!'),end
        k = k-numel(theta(theta==theta(k)));
    end
end
k = 0;
while k<numel(phi)
    k = k+1;
    if phi(k) < 0 || phi(k) > 2*pi
        warning('Out of range phi value ignored!')
        phi = phi(phi~=phi(k));
        if isempty(phi),error('No valid phi value!'),end
        k = k-numel(phi(phi==phi(k)));
    end
end

% generate map for surfacecolor
N = 2^7;
map = jet(N);

npreim = max([numel(theta),numel(phi)]);
surfacecolor = {length(theta)};

% call figure
h = figure(fignum);
if clfIO
    clf
end

% calculate and patch preimages of th/ph point on S2
for preim=1:npreim
    try
        th = theta(preim);
    catch
        th = theta(end);
    end
    try
        ph = phi(preim);
    catch
        ph = phi(end);
    end
    % calculate scalar valued field lim->0 is the preimage of th/ph on S2
    px = sin(th)*cos(ph);
    py = sin(th)*sin(ph);
    pz = cos(th);
    diffmag = sqrt((nn(:,:,:,1)-px).^2+...
                   (nn(:,:,:,2)-py).^2+...
                   (nn(:,:,:,3)-pz).^2);
    % generate an isosurface
    fv = isosurface(diffmag,epsilon);
    % prevents error from an empty patch
    if length(fv.vertices)>=3
        fv.vertices = fv.vertices-...
            repmat([(m-1)/2+1 (n-1)/2+1 (p-1)/2+1],...
            length(fv.vertices),1);
    end
    % determine patch color
    val = round(abs(th)/pi*(N-1)+1);
    surfacecolor{preim} = map(val,:);
    % generate the patch
    hold on
    ptch = patch(fv,...
        'FaceColor',surfacecolor{preim},...
        'FaceAlpha',alpha,...
        'EdgeColor','none',...
        'SpecularStrength',.9,...
        'AmbientStrength',0.3);
    hold off
    ptch.FaceLighting = 'gouraud';
end

axis equal tight off
h.Color = 'w';
view(az,el)
camlight('right')

