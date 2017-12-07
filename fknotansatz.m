% Discription:
% 
% fknotansatz --> fuction that is used to generate a knotted field ansatz.
%   Takes a variable number of arguments in the form of ('flag',input) and
%   outputs a 4D double array with size [nx,ny,nz,3].
%   
%     input: () if nargin = 0, the defalt valus are used
%     optional input: (alpha,beta,a,b,nx,ny,nz,lx,ly,lz)
%       alpha,beta: integer where alpha > beta [2],[1]
%       a,b: integers where the field will form an a-b knot [3],[2]
%       nx,ny,nz: number of points along x, y, and z dimensions [150]
%       lx,ly,lz: length of computational grid [4]
%       f(r): a function f(0)=pi, f(r->inf)=0 ['pi*sech(r)']
%       wknot: change w function for other knots 
%           ['z1.^a./z0.^b'], axis symmetric links [defalt]
%            'z1.^alpha.*z0.^beta./(z1.^a+z0.^b)', trefoil
%     output: (nn)
%       nn: 4D double array with size [nx,ny,nz,3]
% 
%   An example of using this function and visualizing the resulting field:
%       >> clear, nn = fknotansatz; fvisualizePreimages
% 

function nn = fknotansatz(varargin)

% Defalt values:
alpha = 2; beta = 1; % alpha > beta ex: alpha = 2; beta = 1;
a = 3; b = 2; % positive integers ex: a = 3; b = 2;
nx = 150; ny = 150; nz = 150;
lx = 4; ly = 4; lz = 4;
f = @(r) pi*sech(pi*r/2); % f(r) --> f(0)=pi, f(inf)=0 ex: pi*sech(r)
wknot = @(alpha,beta,a,b,z1,z0) z1.^a./z0.^b; % for links

% Argument values
nVarargs = length(varargin);
for k = 1:2:nVarargs
    switch lower(varargin{k}) % case insensitive
        case 'f(r)'
            f = @(r) eval(varargin{k+1});
        case 'wknot'
            wknot = @(alpha,beta,a,b,z1,z0) eval(varargin{k+1});
        case 'alpha'
            alpha = varargin{k+1};
        case 'beta'
            beta = varargin{k+1};
        case 'a'
            a = varargin{k+1};
        case 'b'
            b = varargin{k+1};
        case 'nx'
            nx = varargin{k+1};
        case 'ny'
            ny = varargin{k+1};
        case 'nz'
            nz = varargin{k+1};
        case 'lx'
            lx = varargin{k+1};
        case 'ly'
            ly = varargin{k+1};
        case 'lz'
            lz = varargin{k+1};
    end
end

x1 = linspace(-lx/2,lx/2,nx);
x2 = linspace(-ly/2,ly/2,ny);
x3 = linspace(-lz/2,lz/2,nz);
[x1,x2,x3] = ndgrid(x1,x2,x3);

r = sqrt(x1.^2+x2.^2+x3.^2);
z1 = sin(f(r)).*(x1+1i*x2)./r;
z0 = cos(f(r))+sin(f(r)).*1i.*x3./r;

w = wknot(alpha,beta,a,b,z1,z0);
n1 = (w+conj(w))./(1+abs(w).^2);
n2 = -1i*(w-conj(w))./(1+abs(w).^2);
n3 = (-1+abs(w).^2)./(1+abs(w).^2);

% make nn and force vertical BCs
nn = cat(4,n1,n2,n3);
nn(:,:,1,1) = 0; nn(:,:,1,2) = 0; nn(:,:,1,3) = -1;
nn(:,:,end,1) = 0; nn(:,:,end,2) = 0; nn(:,:,end,3) = -1;
nn(:,1,:,1) = 0; nn(:,1,:,2) = 0; nn(:,1,:,3) = -1;
nn(:,end,:,1) = 0; nn(:,end,:,2) = 0; nn(:,end,:,3) = -1;
nn(1,:,:,1) = 0; nn(1,:,:,2) = 0; nn(1,:,:,3) = -1;
nn(end,:,:,1) = 0; nn(end,:,:,2) = 0; nn(end,:,:,3) = -1;

