function P_AVG = geodesicGroupTest(N,M,option)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N dimension of matrix
% M number of subjects in each group
% option: 1  group action
%         2  parallel transport

%%
P1 = []; P2 = [];
V1 = []; V2 = []; 

%% load all geodesics from the two groups
for k = 1:M
    filename = sprintf('gr1/p_%d.mat',k);
    p = load(filename);
    P1{k} = p.p;
    filename = sprintf('gr2/p_%d.mat',k);
    p = load(filename);
    P2{k} = p.p;
    
    filename = sprintf('gr1/v_%d.mat',k);
    v = load(filename);
    V1{k} = v.v;
    filename = sprintf('gr2/v_%d.mat',k);
    v = load(filename);
    V2{k} = v.v;
end

%% computing the frechet mean
initial_m = zeros(N,N);
for k = 1:M
    filename = sprintf('preprocessed1/sim_subject_%d_1.mat',k);
    covMat = load(filename);
    X1{k} = covMat.covMat;
    initial_m = initial_m + X1{k};
    
    filename = sprintf('preprocessed2/sim_subject_%d_1.mat',k);
    covMat = load(filename);
    X2{k} = covMat.covMat;
    initial_m = initial_m + X2{k};
end

initial_m = initial_m / M / 2;
XALL = [X1,X2];
m = pdMean(XALL,initial_m,M,N);

%% Use a bootstrapped distribution of the mean
sampleNum = 200;
P_SUM = zeros(N,N);

for sidx = 1:sampleNum
    bidx = datasample([1:40],40);
    sampledm = pdMean(XALL(bidx),initial_m,M,N);
    
    % map all tangent vectors to the bootstrapped mean 
    for k = 1:M
        if option == 1 % via Group Action
            VP1{k} = pdGroupAct(V1{k},P1{k},sampledm);
            VP2{k} = pdGroupAct(V2{k},P2{k},sampledm);
        end
        if option == 2 % via Parallel Transport
            VP1{k} = pdParaTrans(V1{k},P1{k},sampledm);
            VP2{k} = pdParaTrans(V2{k},P2{k},sampledm);
        end
    end
    
    % perform two-sample ttests for each matrix entry
    P = zeros(N,N);
    for i = 1:N
        for j = 1:N
            x = zeros(1,M);
            y = zeros(1,M);
            for k = 1:M
                x(k) = VP1{k}(i,j);
                y(k) = VP2{k}(i,j);
            end
        
            [h,P(i,j)] = ttest2(x,y);
        end
    end
    P_SUM = P_SUM + P; 
    P_COL{sidx} = P;
end

P_AVG = P_SUM/sampleNum;

end

