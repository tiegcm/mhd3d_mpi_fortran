%% This script reads .h5 data from mhd3d2

clear all;
close all;

% set(0,'DefaultFigureColormap',jet);

colormap(jet(256));

%% USER INPUTS
% File location
%dirloc = '/global/scratch/bzhu/mhd/otv/pdm_ab_1024/';
dirloc = '/global/scratch/bzhu/test/otv/2048_13d/';
%dirloc = './';
%outdir = '/ihome/bzhu/code/plot/';
%outdir = './';
outdir = dirloc;

% Variables to output
varname = 'p';
dset = strcat('/p');

% Frames saves
nframes = 75;

% dnframes = numel(num2str(nframes)); %number of digits in nframes.
dnframes = 4;
fform = strcat('%0',num2str(dnframes),'d');  % File number format

% Grid size
nx0 = 2048;
ny0 = 2048;
nz0 = 1;

% Box size
alx = 1.;
aly = 1.;
alz = 1.;

%% Stuff derived from user inputs
nx = nx0+8;
ny = ny0+8;

dx = alx/nx0;
dy = aly/ny0;
dz = alz/nz0;

x = [-2:nx0+2]*dx;
y = [-2:ny0+1]*dy;
y0 = [0:ny0-1]*dy;
z = [0:nz0-1]*dz;

ic=25/(36*pi);
file = strcat(dirloc,varname,'.h5');
%dset = strcat('/rho');

hfig=figure(1);

%subplot(2,2,4)

set(hfig,'Position',[1 1 1120 940])
x0=((0:1:nx0-1)+.5)/nx0;
[X,Y]=ndgrid(x0);

fid=fopen(strcat(dirloc,'stime.out'), 'rb');

for n=1:nframes+1
    
    n
    var = h5read(file,dset,[1 1 double(n)],[nx0 ny0 1]);
    p=pcolor(var');
%    p=pcolor(Y,X,var');
    colorbar
    shading interp
%    caxis([.01,.7])

    fseek(fid, 4, 'cof');
    simt = fread(fid,[1,1],'float64');
    fseek(fid, 4, 'cof');
    
    title(strcat('t=',num2str(simt)))
%    title(strcat('n=',num2str(n-1)))

    
%    plot(x0,var(:,1),'-o');
%    plot(var(1,:),'-o');

%    xlim([.9,1.])
%    print(hfig,'-dpng',strcat(outdir,varname,num2str(n-1+n0,fform)));

    pause(.2)

end

%% test symmetry
for i=1:nx0
    for j=1:ny0
        diff(i,j)=var(i,j)-var(nx0+1-i,ny0+1-j);
    end
end
figure(2)
p=pcolor(diff');

%fclose(fid);
