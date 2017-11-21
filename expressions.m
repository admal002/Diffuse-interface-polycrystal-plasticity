% Copyright (c) 2017 Nikhil Admal/Jaime Marian
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy of
% this software and associated documentation files (the 'Software'), to deal in
% the Software without restriction, including without limitation the rights to
% use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
% the Software, and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
  
% THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
% 
% Author: Nikhil Admal
% 
% 
% This matlab script is used to generate two .txt files: 
% 1) variables.txt and 2) pde.txt
% The above files are subsequently used to run a comsol implementation of 
% a three-dimensional simulation of a diffuse-interface polycrystality 
% model as described in the paper:
%
%   N. C. Admal, G. Po, J. Marian. A unified framework for polycrystal 
%   plasticity with grain boundary evolution.
%   arXiv preprint arXiv:1709.10176, 2017. 
%   Submitted to International Journal of Plasticity

% Additional details on the usage of the above three files for the
% implementation of COMSOL are given in the README file. 

clear all
close all
clc

dim=3;
numSlips=48;


fID=fopen('variables.txt','w');
lineFormat='%s\n';

%% u
u={'u1','u2','u3'};
X={'x','y','z'};

%% w_el
syms c11;
syms c12;
syms c44;

%% Define symbolic components of Ee, Ue
for bi=1:dim
    for bj=1:dim
        eval(['syms Ee' num2str(bi) num2str(bj) ' real']);
    end
end

%% Assemble the symbolic matrix Ee

Ee = [Ee11, Ee12, Ee13; 
      Ee21, Ee22, Ee23;
      Ee31, Ee32, Ee33];
 
%% Energy

% w_el is a function of Ee

Cijkl=[c11, c12, c12,   0,   0,   0;
       c12, c11, c12,   0,   0,   0;
       c12, c12, c11,   0,   0,   0;
       0,   0,   0,   c44,   0,   0;
       0,   0,   0,     0, c44,   0;
       0,   0,   0,     0,   0, c44];

ee=[Ee11, Ee22, Ee33, Ee23+Ee32, Ee13+Ee31, Ee12+Ee21];

w=0;
for i=1:6
    for j=1:6
        w=w+Cijkl(i,j)*ee(j)*ee(i);
    end
end
w=0.5*w;


varName='w';
description=['Elastic energy'];
expression = char(eval(w));
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);



%% Derivatives of energy
%  dw/dEe


% W (dw/dEe)
% W_{bi,bj} = d(w,Ee_{bi,bj})
for bi=1:dim
    for bj=1:dim
        varName=['W' num2str(bi) num2str(bj)];
        description=['intermediate stress, component ' num2str(bi) num2str(bj)];
        expression=[];
        expression = char(eval(['diff(w,Ee' num2str(bi) num2str(bj) ')']));
        expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end



%% Rp

syms q1 real;
syms q2 real
syms q3 real;
syms qt1 real;
syms qt2 real;
syms qt3 real;

psi=[q1,q2,q3];
psiNorm = sqrt(q1^2+q2^2+q3^2);

W=sym(zeros(dim,dim));
for i=1:dim
    for j=1:dim
        for k=1:dim
         W(i,j) = W(i,j)-levCiv(i,j,k)*psi(k);
        end
    end
end

Rp = eye(dim)+(sin(psiNorm)/psiNorm)*W+0.5*(sin(0.5*psiNorm)/(0.5*psiNorm))^2*W*W

dotRp = diff(Rp,q1)*qt1+diff(Rp,q2)*qt2+diff(Rp,q3)*qt3;
  
for bi=1:dim
    for I=1:dim
        varName=['Rp' num2str(bi) num2str(I)];
        description=['Rp, component ' num2str(bi) num2str(I)];
        expression = [char(Rp(bi,I))];
        expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end

%% Up
syms Up11;
syms Up12;
syms Up13;
syms Up22;
syms Up23;
syms Up33;
syms Upt11;
syms Upt12;
syms Upt13;
syms Upt22;
syms Upt23;
syms Upt33;

Up = [Up11, Up12, Up13;
      Up12, Up22, Up23;
      Up13, Up23, Up33];
dotUp = [Upt11, Upt12, Upt13;
         Upt12, Upt22, Upt23;
         Upt13, Upt23, Upt33];



%% Fp = Rp Up
Fp = Rp*Up;
for bi=1:dim
    for I=1:dim
        varName=['Fp' num2str(bi) num2str(I)];
        description=['Fp, component ' num2str(bi) num2str(I)];
        expression = char(Fp(bi,I));
        expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end

%% dot Fp = dotRp*Up + Rp*dotUp

dotFp = dotRp*Up+Rp*dotUp;


%%%%%%%%%



%% da
% The l.h.s of the flow rule \dot \Fp = Lp Fp is implemented as the matrix da 
% times the time-derivative of the unknown array 
% [q1, q2, q3, Up11, Up22, Up33, Up12, Up13, Up23]^T
da=sym(zeros(dim^2,dim^2));

for i=1:dim
    for j=1:dim
        k=(i-1)*dim+j;
        da(k,1)=diff(dotFp(i,j),qt1);
        da(k,2)=diff(dotFp(i,j),qt2);
        da(k,3)=diff(dotFp(i,j),qt3);
        da(k,4)=diff(dotFp(i,j),Upt11);
        da(k,5)=diff(dotFp(i,j),Upt22);
        da(k,6)=diff(dotFp(i,j),Upt33);
        da(k,7)=diff(dotFp(i,j),Upt12);
        da(k,8)=diff(dotFp(i,j),Upt13);
        da(k,9)=diff(dotFp(i,j),Upt23);
    end
end


for i=1:dim*dim
    for j=1:dim*dim
        varName=['da' num2str(i) num2str(j)];
        description=['da, component ' num2str(i) num2str(j)];
        expression=char(da(i,j));
        expression=expression(setdiff([1:length(expression)],findstr(expression,' ')));
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end
%%%%%%%%%%%%%%%%%


%% Gp
for bi=1:dim
    for I=1:dim
        eval(['syms Fp' num2str(bi) num2str(I) ' real;']);
        eval(['Fps(' num2str(bi) ',' num2str(I) ')= Fp' num2str(bi) num2str(I) ';']);
    end
end
Jp=det(Fps);
expression=char(Jp);
expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
fprintf(fID,lineFormat,['Jp' ' ' expression ' ' 'determinant of Fp']);

Gp=inv(Fps);
for bi=1:dim
    for I=1:dim
        varName=['Gp' num2str(I) num2str(bi)];
        description=['inverse Fp, component ' num2str(I) num2str(bi)];
        expression=char(Gp(I,bi)*Jp);
        expression=expression(setdiff([1:length(expression)],strfind(expression,' ')));
        fprintf(fID,lineFormat,[varName ' (' expression ')/Jp ' description]);
    end
end

%% F
% F_{i,J} = d(u_i,X{J}) + delta_{i,J}
for i=1:dim
    for J=1:dim
        varName=['F' num2str(i) num2str(J)];
        description=['deformation gradient, component ' num2str(i) num2str(J)];
        expression=[u{i} X{J} '+' num2str(i==J)];
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end


%% Fe = F*Gp
% Fe_{i,bi} = F_{i,J}*Gp_{J,bi}
for i=1:dim
    for bi=1:dim
        varName=['Fe' num2str(i) num2str(bi)];
        description=['elastic distortion, component ' num2str(i) num2str(bi)];
        expression=[];
        for J=1:dim
            expression=[expression '+F' num2str(i) num2str(J) '*Gp' num2str(J) num2str(bi)];
        end
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end


%% Ce_{bi,bj}=Fe_{k,bi} Fe_{k,bj}
for bi=1:dim
    for bj=1:dim
        varName=['Ce' num2str(bi) num2str(bj)];
        description=['left CG stretch tensor, component ' num2str(bi) num2str(bj)];
        expression=[];
        for k=1:dim
            expression=[expression '+Fe' num2str(k) num2str(bi) '*Fe' num2str(k) num2str(bj)];
        end
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end

%% Ee_{bi,bj} = 0.5*Ce_{bi,bj} - 0.5 delta_{bi,bj}

for bi=1:dim
    for bj=1:dim
        varName=['Ee' num2str(bi) num2str(bj)];
        description=['elastic strain, component ' num2str(bi) num2str(bj)];
        expression= ['0.5*Ce' num2str(bi) num2str(bj) '-0.5*' num2str(bi==bj)];
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end



%% Piola P
% P_{k,I}=Fe_{k,bi}*W_{bi,bj}*Gp_{I,bj}
for k=1:dim
    for I=1:dim
        varName=['P' num2str(k) num2str(I)];
        description=['first Piola stress, component ' num2str(k) num2str(I)];
        expression=[];
        for bi=1:dim
            for bj=1:dim
                expression=[expression '+Fe' num2str(k) num2str(bi) '*W' num2str(bi) num2str(bj) '*Gp' num2str(I) num2str(bj)];
            end
        end
        
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end



%% Resolved shear stress
% rss = W_{bi,bj}*n_{bj}*Ce_{bi,bk}*s_{bk} +
%       \
for slipID=1:numSlips
    varName=['rss' num2str(slipID)];
    description=['Resolved shear stress ' num2str(slipID)];
    expression=[];
    % W_{bi,bj}*n_{bj}*Ce_{bi,bk}*s_{bk}
    for bi=1:dim
        for bj=1:dim
            for bk=1:dim
                expression=[expression '+W' num2str(bi) num2str(bj) '*n' num2str(slipID) '_' num2str(bj) '*Ce' num2str(bi) num2str(bk) '*s' num2str(slipID) '_' num2str(bk)];
            end
        end
    end
    fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
end

%% vi = v0*(taui/g)^m
for slipID=1:numSlips
    varName=['v' num2str(slipID)];
    description=['slip rate'];
    expression=['sign(rss' num2str(slipID) ')*v0*(abs(rss' num2str(slipID) ')/g)^m'];
    fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
end

%% g = g0 + b*c44*sqrt(hn*rhon)
varName='g';
description=['yield stress'];
expression='g0+b*c44*sqrt(hn*rhon)';
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

%% rhon = rho0*h
varName='rhon';
description=['dislocation density'];
expression='rho0*h';
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

%% v
varName='v';
description=['total slip rate'];
expression=[];
for slipID=1:numSlips
    expression=[expression '+abs(v' num2str(slipID) ')'];
    %expression=[expression '+v' num2str(slipID) '^2'];
end
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

%% doth = v*k1*sqrt(h)-k2*h
varName='doth';
description=['rate of change in dislocation density'];
expression='v*k1*sqrt(h)-k2*h';
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);

%% k2 = k20*(vk0*v^(n-1))(1/n)
varName='k2';
description=['dislocation annihilation'];
expression='k20*vk0^(1/n)*v^((n-1)/n)';
fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);


%% Lp
for bi=1:dim
    for bj=1:dim
        varName=['Lp' num2str(bi) num2str(bj)];
        description=['Lp, component ' num2str(bi) num2str(bj)];
        expression=[];
        for slipID=1:numSlips
            expression=[expression '+v' num2str(slipID) '*s' num2str(slipID) '_' num2str(bi) '*n' num2str(slipID) '_' num2str(bj)];
        end
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end


%% LpFp
for bi=1:dim
    for I=1:dim
        varName=['LpFp' num2str(bi) num2str(I)];
        description=['LpFp, component ' num2str(bi) num2str(I)];
        expression=[];
        for bj=1:dim
            expression=[expression '+Lp' num2str(bi) num2str(bj) '*Fp' num2str(bj) num2str(I)];
        end
        fprintf(fID,lineFormat,[varName ' ' expression ' ' description]);
    end
end

%% G 
% G = Fp_{bi,K} * eps(K,R,S)*d(Fp_{bj,S},R)
for bi=1:3
    for bj=1:3
        varName=['G' num2str(bi) num2str(bj)];
        description=['G, component ' num2str(bi) num2str(bj)];
        expression=[];
        for K=1:3
            for R=1:2
                for S=1:3
                    expression=[expression '+' num2str(levCiv(K,R,S)) '*d(Fp' num2str(bj) num2str(S) ',' X{R} ')*Fp' num2str(bi) num2str(K)];
                end
            end
        end
            fprintf(fID,lineFormat,[varName ' (' expression ')/Jp ' description]);
    end
end


%% Close the file

fclose(fID);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PDE for crystal plasticity

FID=fopen('pde.txt','w');


% P_{iJ} * test(ui,X{J})
for i=1:dim
    expression=[];
    for J=1:dim
        expression=[expression '+P' num2str(i) num2str(J) '*test(u' num2str(i) X{J} ')'];
    end
    fprintf(FID,[expression '\n']);
end

fclose(fID);


