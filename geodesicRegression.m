function geodesicRegression(M)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each subject, load 3 subject-specific matrices
% compute geodesic regression to get a starting point and a shooting
% direction


%% group 1
for k = 1:M
    filename1 = sprintf('preprocessed1/sim_subject_%d_1.mat',k);
    covMat = load(filename1);
    X{1} = covMat.covMat;
    filename2 = sprintf('preprocessed1/sim_subject_%d_2.mat',k);
    covMat = load(filename2);
    X{2} = covMat.covMat;
    filename3 = sprintf('preprocessed1/sim_subject_%d_3.mat',k);
    covMat = load(filename3);
    X{3} = covMat.covMat;
    
    [p,v] = pdRegression(X);
    
    filename = sprintf('gr1/p_%d.mat',k);
    save(filename,'p');
    filename = sprintf('gr1/v_%d.mat',k);
    save(filename,'v');
    
    fprintf('**Finish Geodesic Regression Subject %d ctrl\n',k);
end


%% group 2
for k = 1:M
    filename1 = sprintf('preprocessed2/sim_subject_%d_1.mat',k);
    covMat = load(filename1);
    X{1} = covMat.covMat;
    filename2 = sprintf('preprocessed2/sim_subject_%d_2.mat',k);
    covMat = load(filename2);
    X{2} = covMat.covMat;
    filename3 = sprintf('preprocessed2/sim_subject_%d_3.mat',k);
    covMat = load(filename3);
    X{3} = covMat.covMat;
    
    [p,v] = pdRegression(X);
    
    filename = sprintf('gr2/p_%d.mat',k);
    save(filename,'p');
    filename = sprintf('gr2/v_%d.mat',k);
    save(filename,'v');
    
    fprintf('**Finish Geodesic Regression Subject %d disease\n',k);
end

end