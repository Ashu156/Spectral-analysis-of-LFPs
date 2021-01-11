%%

%% %% Prepare dataset

% [fileName, path] = uigetfile;
% load(strcat(path, fileName));

labels = cell(size(data, 1), 1);

for i = 1:size(data, 1)
    
    if i <= size(control_total_data_frex_chunkwise, 1)

        labels{i} = 'control';
    else
        labels{i, 1} = 'stress';
    end
end

labels_num = grp2idx(labels); % change labels to numeric indices

%%

% binary classification 
% X = randn(100,10); % Adding some random features
% X(:,[1,3,5,7]) = meas(1:100,:); % 1, 3, 5, 7 are actula features in the dataset while 2, 4, 6, 8, 9, 10 ae dummy faetures
% y = labels_num(1:100); % labels for the instances

%% Partitioning the training and the testing dataset (80-20 partition here)

X = data;
y = labels_num;
train_proportion = 0.8;
nIter = 1e3;
test_accuracy = zeros(10, 1);

 for ii = 1:nIter
     
     rand_num = randperm(size(X,1)); % shuffled integers
     X_train = X(rand_num(1:round(train_proportion*length(rand_num))),:); % training data
     y_train = y(rand_num(1:round(train_proportion*length(rand_num))),:); % training labels

      X_test = X(rand_num(round(train_proportion*length(rand_num))+1:end),:); % testing data
      y_test = y(rand_num(round(train_proportion*length(rand_num))+1:end),:); % testing labels

% CV partition

      c = cvpartition(y_train,'k', 10); % partioning training data for 5-fold cross-validation

      Md1 = fitcsvm(X_train,y_train,'KernelFunction', 'linear', 'OptimizeHyperparameters', 'auto',...
      'HyperparameterOptimizationOptions', struct('AcquisitionFunctionName',...
      'expected-improvement-plus','ShowPlots', false)); % Bayes' Optimization

       test_accuracy_for_iter = sum((predict(Md1, X_test) == y_test))/length(y_test)*100; % accuracy
       test_accuracy(ii, 1) = test_accuracy_for_iter;
       disp(strcat('Iteration_ ', num2str(ii), '_running'))
 end

%% hyperplane

% figure;
% hgscatter = gscatter(X_train_w_best_feature(:,1),X_train_w_best_feature(:,2),y_train);
% hold on;
% h_sv = plot(Md1.SupportVectors(:,1),Md1.SupportVectors(:,2),'ko','markersize',8);
% 
% 
% % test set? data? ?? ??? ????.
% 
% gscatter(X_test_w_best_feature(:,1),X_test_w_best_feature(:,2),y_test,'rb','xx')
% 
% % decision plane
% XLIMs = get(gca,'xlim');
% YLIMs = get(gca,'ylim');
% [xi,yi] = meshgrid([XLIMs(1):0.01:XLIMs(2)],[YLIMs(1):0.01:YLIMs(2)]);
% dd = [xi(:), yi(:)];
% pred_mesh = predict(Md1, dd);
% redcolor = [1, 0.8, 0.8];
% bluecolor = [0.8, 0.8, 1];
% pos = find(pred_mesh == 1);
% h1 = plot(dd(pos,1), dd(pos,2),'s','color',redcolor,'Markersize',5,'MarkerEdgeColor',redcolor,'MarkerFaceColor',redcolor);
% pos = find(pred_mesh == 2);
% h2 = plot(dd(pos,1), dd(pos,2),'s','color',bluecolor,'Markersize',5,'MarkerEdgeColor',bluecolor,'MarkerFaceColor',bluecolor);
% uistack(h1,'bottom');
% uistack(h2,'bottom');
% legend([hgscatter;h_sv],{'control','stress','support vectors'})

%% end of script