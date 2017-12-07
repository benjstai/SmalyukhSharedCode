% Discription:
% 
% fvisualizePreimages --> a fuction that is used to visualize preimages, it
%   takes a variable number of arguments in the form of -flag [argument].
%   Number of preimages plotted is max([numel(theta),numel(phi)]).
%     If theta(k) or phi(k) exceeds the index, the last theta/phi valus is
%     used eg. theta(end) or phi(end).
%   Lighing is applied to each surface during creation.
%   Use ('clf',false) to use the figure without clearing (defalt is true).
% 
%     input: () if nargin = 0, nn from the base workspace is used
%     optional input:  (nn, theta, phi, alpha, epsilon, fignum, clf)
%       nn: 4D double array with size [m,n,p,3]
%       theta: list of polar angles to use
%       phi: list of azamuthal angle to use
%       alpha: transparency of preimages
%       epsilon: value used to generate isosurface in 'diffmag'
%       fignum: figure to use
%       clf: toggle clf true/false
% 

function fvisualizePreimages(varargin)

% Defalt values:
nn = [];
theta = 0;
phi = 0;
alpha = 1;
epsilon = 0.3;
fignum = 1;
clfIO = 1;

% new defalt value
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
    end
end

% if nn is not given, look for nn in the base workspace
if nargin==0 || isempty(nn)
    try
        disp('Looking for ''nn'' in the base Workspace')
        nn = evalin('base','nn');
    catch
        error('no ''nn'' variable found in the base Workspace')
    end
end

[m,n,p,test] = size(nn);
if ~(test==3)
    error('First argument must be a 4D double array with size [m,n,p,3]')
end

% if numel(theta)~=numel(phi)
%     error('numel(theta)~=numel(phi)')
% end

% generate map for surfacecolor
N = 2^7;
map = jet(N);

npreim = max([numel(theta),numel(phi)]);
surfacecolor = {length(theta)};
h = figure(fignum);
if clfIO
    clf
end
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
    px = sin(th)*cos(ph);
    py = sin(th)*sin(ph);
    pz = cos(th);
    diffmag = sqrt((nn(:,:,:,1)-px).^2+...
                   (nn(:,:,:,2)-py).^2+...
                   (nn(:,:,:,3)-pz).^2);
    fv = isosurface(diffmag,epsilon);
    if length(fv.vertices)>=3
        fv.vertices = fv.vertices-...
            repmat([(m-1)/2+1 (n-1)/2+1 (p-1)/2+1],...
            length(fv.vertices),1);
    end
    val = round(abs(th)/pi*(N-1)+1);
    surfacecolor{preim} = map(val,:);
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
camlight('right')
