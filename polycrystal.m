%
%%%%%%%%%%%%%%%%%%%%%%%%%%   Polycystal generator  %%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Nikhil Admal                                                      %
% Institution: University of California Los Angeles                         %
% Department: Materials Science and Engineering                             %
%
% A polycrystal is viewed as a Voronoi tessellation associated with a
% distribution of grain centers and their corresponding random
% orientations.
%
% The program constructs grain centers, and their corresponding random
% orientations for a given 
% log-normal distribution for grain sizes. Size of a grain is defined 
% as the minimum of the distances between the its center to
% other grain enters.
%
% Program output: The final polycrystal is not outputted as a Voronoi tessellation.
% Instead, it is outputted in the form of a uniform grid, with each grid
% point associated to the closes grid center. Additionally the grains may be 
% scaled to form plates/needle-like microstructures.


clear all
close all
clc

% The code has not been tested for dim = 2.
% So do not change the value of dim
dim = 3;

% Scales used to generate different textures such as needles or plates
scalex = 1.0;
scaley = 1.0;
scalez = 1.0;


% nGrains: Number of grains
% Note: This does not mean we end up with 'nGrains' grains.
%       The program tries to pack 'nGrains' grains subjected to a given 
%       distribution of grain size over 'maxit' iterations.

nGrains = 150;

% length of the sample
length = 3e-6;

% Maximum twist about an axis of rotation
maxtwist = 20*pi/180;

% Covergence 
maxit = 1000;

% Resolution of the grid to output the polycrystal
ngrid=100;

% The algorithm begins here.

% Variables
% ---------
%     ctr...................Counter for number of grains
%     grainSize.............An array containing the grain size normalized w.r.t to length
%     randPoint.............A random point in the domain chosen from a uniform distribution
%     q.....................A random vector describing the orientation of the grain.
%                           |q| = twist, and q/|q| = dir which is the rotation axis. 
%                           dir is sampled from a uniform distribution on a sphere.
%     grain(dimn,nGrans)....Grain centers
%     b[]...................Array with entries distribute log normally
%         
   

ctr = 1;
grainSize = zeros(nGrains,1);

b=[];
% Iterate until maxit
for i=1:maxit
    % Construct a random point from a uniform distibuion 
    for j=1:dim
        randPoint(j) = rand;
    end
    
    % Construct a random orientation vector q
    ru=rand;
    theta=rand*pi/2;
    dir = [sqrt(1-ru^2)*cos(theta);...
           sqrt(1-ru^2)*sin(theta);...
           ru];
    twist = rand*maxtwist;

    q=twist*dir;
    
   % If it the first random point, then accept it as the first grain
    if (ctr == 1) 
        grain(:,ctr) = randPoint;
        
        orientation(ctr,:) = q;
        if (ctr == nGrains)
                break;
        end
        ctr = ctr+1;
    else
        % If it is not the first grain center, then measure its distrance to
        % existing grain centers and store in the array dist
        dist=zeros(ctr-1,1);

        for k=1:ctr-1
            dist(k)=(randPoint(1)-grain(1,k))^2/scalex^2+...
                    (randPoint(2)-grain(2,k))^2/scaley^2+...
                    (randPoint(3)-grain(3,k))^2/scalez^2;
        end
        dist = sqrt(dist);
        
        % rpoint is a random variable (rv) with a uniform distribution, used to
        % generate a rv logcdf^{-1}(x), with a log-normal distribution,
        % where logcdf is the cdf of a log-normal distribution
        % Note: Keep rpoint away from 0 and 1 as it gets difficult to find
        % the roots of logcdf^{-1}(x)-rpoint=0, when rpoints is close to 0
        % or 1
        rpoint=(0.995-0.0015)*rand+0.0015;


        fun=@(x) logcdf(x)-rpoint;
        
        % minDist is a random variable with log-normal distribution
        minDist=fzero(fun,0.2);
        b=[b;minDist];
        
        % Accept randPoint as the grain center if all the distances to
        % existing grid centers are greater than minDist.
        if (all(dist > minDist))
            grain(:,ctr) = randPoint;
            
            % Once randPoint is accepted as a grid center, associate the
            % randorm orientation to it.
            orientation(ctr,:) = q;
            
            % If we have reached the 'nGrains' grains, then exit. Else
            % iterate.
            if (ctr == nGrains)
                break;
            end
            
            ctr = ctr + 1
        end
    end
    if (i == maxit)
        nGrains = ctr-1;
    end
end

% Calculate the grain sizes
for i = 1:nGrains
    dist=[];
    for j = 1:nGrains
        if (j ~= i)
                dist = [dist; (grain(1,i)-grain(1,j))^2+...
                              (grain(2,i)-grain(2,j))^2+...
                              (grain(3,i)-grain(3,j))^2];
        end
    end
    grainSize(i) = sqrt(min(dist));
end
            
% Compare the pdfs of the actual log normal distribution using 'b' and the
% distribution of grain sizes using 'grainSize' in a histogram plot
figure;
histogram(grainSize(1:nGrains),'BinWidth',0.01,'Normalization','pdf')
hold on;
histogram(b,'BinWidth',0.01,'Normalization','pdf')

tix=get(gca,'xtick')';
set(gca,'xticklabel',num2str(tix,'%.2f'))
tix=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tix,'%.2f'))

% Ourput the voronoi tessellation in the form of a grid
x=[0:1/ngrid:1];
y=[0:1/ngrid:1];
if (dim==2)
    [X,Y] = meshgrid(x,y);
elseif (dim==3)
    z=[0:1/ngrid:1];
    [X,Y,Z] = meshgrid(x,y,z);
end

maxGrainPoints = ngrid^3;
numGrainPoints = zeros(nGrains,1);
for i=1:numel(X)
    if (dim == 2) 
        point = [X(i) Y(i)];
    elseif (dim==3)        
        point = [X(i) Y(i) Z(i)];
    end
    dist = zeros(nGrains,1);
    for j=1:nGrains
        for k=1:dim
            dist(j) = (point(1)-grain(1,j))^2/scalex^2+...
                      (point(2)-grain(2,j))^2/scaley^2+...
                      (point(3)-grain(3,j))^2/scalez^2;
        end
    end
    dist=sqrt(dist);
    [minVal,minI] = min(dist);
    F1(i) = orientation(minI,1);
    F2(i) = orientation(minI,2);
    F3(i) = orientation(minI,3);
     numGrainPoints(minI) = numGrainPoints(minI)+1;
    if (numGrainPoints(minI) > maxGrainPoints)
        numGrainPoints(minI)
        break;
    end

end

F1 = reshape(F1,size(X));
F2 = reshape(F2,size(X));
F3 = reshape(F3,size(X));

% Output the polycrystal in terms of the grid point and associated
% orientation. The output is taylored to be an input file for COMSOL

fid1 = fopen('orientation1.txt','w');
fid2 = fopen('orientation2.txt','w');
fid3 = fopen('orientation3.txt','w');

gridArray = length*[0:1/ngrid:1];
fprintf(fid1,'%s\n', '%grid');
fprintf(fid1,'%.5e ', gridArray);
fprintf(fid1,'\n');
fprintf(fid1,'%.5e ', gridArray);
fprintf(fid1,'\n');
fprintf(fid1,'%.5e ', gridArray);
fprintf(fid1,'\n');
fprintf(fid1,'%s', '%data');

fprintf(fid2,'%s\n', '%grid');
fprintf(fid2,'%.5e ', gridArray);
fprintf(fid2,'\n');
fprintf(fid2,'%.5e ', gridArray);
fprintf(fid2,'\n');
fprintf(fid2,'%.5e ', gridArray);
fprintf(fid2,'\n');
fprintf(fid2,'%s', '%data');

fprintf(fid3,'%s\n', '%grid');
fprintf(fid3,'%.5e ', gridArray);
fprintf(fid3,'\n');
fprintf(fid3,'%.5e ', gridArray);
fprintf(fid3,'\n');
fprintf(fid3,'%.5e ', gridArray);
fprintf(fid3,'\n');
fprintf(fid3,'%s', '%data');



for i=1:ngrid+1
    fprintf(fid1,'\n');
    fprintf(fid2,'\n');
    fprintf(fid3,'\n');
    for j=1:ngrid+1
        for k=1:ngrid+1
            fprintf(fid1,'%6.4f ',F1(k,i,j));
            fprintf(fid2,'%6.4f ',F2(k,i,j));
            fprintf(fid3,'%6.4f ',F3(k,i,j));
        end
    end
end
fclose(fid1);
fclose(fid2);
fclose(fid3);

% figure
% 
% for i=1:3
%     
%     scatter3(grainPoints(i,:,1),grainPoints(i,:,2),grainPoints(i,:,3),50,orientation(i,1)*ones(1,maxGrainPoints),'filled');
%     xlim([0 1]);
%     ylim([0 1]);
%     zlim([0 1]);
% 
%     hold on;
% end
% 
% scatter3(grain(1,1:3),grain(2,1:3),grain(3,1:3),100,'red','filled');
% xlim([0 1]);
% ylim([0 1]);
% zlim([0 1]);
% hold on;
% 
% 
% 
% % Plot 
% % figure
% % scatter3(pointsArray(1,:),pointsArray(2,:),pointsArray(3,:),25,orientation(:,1),'filled');
% % hold on;
% 
% bx = reshape(X(:,:,1),[],1);
% by = reshape(Y(:,:,1),[],1);
% bz = reshape(Z(:,:,1),[],1);
% bF = reshape(F1(:,:,1),[],1);
% 
% tx = reshape(X(:,:,ngrid+1),[],1);
% ty = reshape(Y(:,:,ngrid+1),[],1);
% tz = reshape(Z(:,:,ngrid+1),[],1);
% tF = reshape(F1(:,:,ngrid+1),[],1);
% 
% lx = reshape(X(1,:,:),[],1);
% ly = reshape(Y(1,:,:),[],1);
% lz = reshape(Z(1,:,:),[],1);
% lF = reshape(F1(1,:,:),[],1);
% 
% rx = reshape(X(ngrid+1,:,:),[],1);
% ry = reshape(Y(ngrid+1,:,:),[],1);
% rz = reshape(Z(ngrid+1,:,:),[],1);
% rF = reshape(F1(ngrid+1,:,:),[],1);
% 
% fx = reshape(X(:,1,:),[],1);
% fy = reshape(Y(:,1,:),[],1);
% fz = reshape(Z(:,1,:),[],1);
% fF = reshape(F1(:,1,:),[],1);
% 
% dx = reshape(X(:,ngrid+1,:),[],1);
% dy = reshape(Y(:,ngrid+1,:),[],1);
% dz = reshape(Z(:,ngrid+1,:),[],1);
% dF = reshape(F1(:,ngrid+1,:),[],1);
% 
% scatter3(bx,by,bz,50,bF,'filled');
% hold on;
% scatter3(tx,ty,tz,50,tF,'filled');
% hold on;
% scatter3(lx,ly,lz,50,lF,'filled');
% hold on;
% scatter3(rx,ry,rz,50,rF,'filled');
% hold on;
% scatter3(fx,fy,fz,50,fF,'filled');
% hold on;
% scatter3(dx,dy,dz,50,dF,'filled');

        
        