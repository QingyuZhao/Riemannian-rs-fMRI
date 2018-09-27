clear;
clf;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A Riemannian Framework for Longitudinal Analysis of Resting-State
% Functional Connectivity
%
% * A synthetic example of performing elementwise two-sample t-tests on the
% covariance matrices of two groups of subjects;
% * Each subject has 3 covariance matrices, from which the longitudinal
% trajectory is estimated;
% * Detailed experimental explaination given in the following paper
%
% Citation
% A Riemannian Framework for Longitudinal Analysis of Resting-State
% Functional Connectivity, Q Zhao, D Kwon, KM Pohl, International
% Conference on Medical Image Computing and Computer-Assisted Intervention
% 2018
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N dimension of matrix
% M number of subjects in each group
% option: 1  group action
%         2  parallel transport

M = 20;
N = 10;
option = 1;

% perform geodesic regression for each subject
geodesicRegression(M);
% perform group analysis on the resulting tangent vectors
p = geodesicGroupTest(N,M,option);
% visualize p-value map
imagesc(p);
title('p-value map');
colorbar;

