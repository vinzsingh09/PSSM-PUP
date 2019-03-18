function foldresult = ComparePerformanceAssessment( testdata, comparedata)

	testlabelY = cell2mat(testdata(:,2));
    %comparedatalabelY = cell2mat(comparedata(:,2));
    comparedatalabelY = comparedata;
	confumat = confusionmat(testlabelY, comparedatalabelY);

	TN = confumat(1,1);
	FP = confumat(1,2);
	FN = confumat(2,1);
	TP = confumat(2,2);

	Sensitivity = (TP / (TP + FN)) * 100;
	Specificity = (TN / (TN + FP)) * 100;
	Precision = (TP / (TP + FP)) * 100;
	Accuracy = ((TP + TN)/ (TP + FN + FP + TN))*100;
	MCC = ((TN * TP) - (FN * FP)) / sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN));


	%folds	TP	FN	FP	TN	Sensitivity (Sn)	Specificity (Sp)	Precision (Pre)	Accuracy (ACC)	MCC	AUC

	foldresult = [TP FN	FP	TN	Sensitivity Specificity	Precision Accuracy MCC];
end

