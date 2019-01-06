%------------------------------ RQNN-PCA ----------------------------------
% Author=Bruno Losseau                                      Date=6/01/2019
% Implementation of Pompili's never published algorithm P-ONPMF-R using
% Manopt and with the canvas of the patch 
% Optimization-on-manifolds-with-extra-constraints-master by Liu and Boumal
% available at https://github.com/losangle.
%
% only Manopt must be on your computer and linked to matlab using the
% setpath (available at https://www.manopt.org/)

%close all; clear all; clc;
specifier.matlabversion = 1; %0 if older than 2015 1 otherwise


% We provide an example with histology images. M is a RGB image containing
% nuclei. We proceed to its segmentation based on 2 clusters
% images source: https://www.dropbox.com/s/9knzkp9g9xt6ipb/cluster%20nuclei.zip?dl=0


load('histologyB.mat'); 

%_______Set up data______
k=numClasses;
M = M1;
M(M<0)=0;

%________Experiment_____
options.maxOuterIter = 2000;
options.maxtime = 3600;
options.minstepsize = 1e-5  ;
options.verbosity=0;
options.showErr=0;
tic()
[result,U] = clientconstraint_stiefel_NNPCA(M, k, options);
toc()

figure()

showImage(U,k,[sz,sz],1,1)
