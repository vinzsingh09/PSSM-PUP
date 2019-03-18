function foldresult = PerformanceAssessment( traindata, testdata)
    dataX=cell2mat(traindata(:,1));
	labelY = cell2mat(traindata(:,2));
	%nTrees = 50;
	%for i=1:length(traindata(:,2))
	%labelY{i,1} = int2num(labelY{i});
	%end
	%PrediciveModel = TreeBagger(nTrees,dataX,lab elY, 'Method', 'classification');
    weight_vector = randn (length(labelY), 1);
    %model = svmtrain([], labelY, dataX, '-t 2 -c 1 -g 0.050 -r 3');
    %model = svmtrain([], labelY, dataX, '-t 1 -c 1 -g 0.025 -r 3');
     model = svmtrain([], labelY, dataX, '-t 2 -c 2 -g 0.0250 -b 1 -w1 1 -w-1 1');
	
%oobPredict(PrediciveModel);
	testX = cell2mat(testdata(:,1));
	%predChar = PrediciveModel.predict(testX);
	testlabelY = cell2mat(testdata(:,2));
    [predict_label, accuracy, dec_values] = svmpredict(testlabelY, testX, model);

	%for i=1:length(testdata(:,2))
	%testlabelY{i,1} = int2str(testlabelY{i});
	%end

	confumat = confusionmat(testlabelY, predict_label);

	TN = confumat(1,1);
	FP = confumat(1,2);
	FN = confumat(2,1);
	TP = confumat(2,2);

	Sensitivity = (TP / (TP + FN)) * 100;
	Specificity = (TN / (TN + FP)) * 100;
	Precision = (TP / (TP + FP)) * 100;
	Accuracy = ((TP + TN)/ (TP + FN + FP + TN))*100;
	MCC = ((TN * TP) - (FN * FP)) / sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN));
    %gmean = sqrt(Sensitivity * Specificity);


	%folds	TP	FN	FP	TN	Sensitivity (Sn)	Specificity (Sp)	Precision (Pre)	Accuracy (ACC)	MCC	AUC

	foldresult = [TP FN	FP TN Sensitivity Specificity Precision Accuracy MCC];
    %foldresult = [TP FN	FP	TN	 Accuracy ];
end

