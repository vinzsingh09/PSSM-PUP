clear 

load 'C:\Users\singh_vn\Desktop\MSc\Data\Pupylation\Matlab Data\Pupylationstruct.mat'

sequence_size = length(DB_Pupylation);

% counting K  // 2471 K
% count_k = 2471;
count_k = 0;
for i=1:sequence_size
    count_k = count_k + length(strfind(DB_Pupylation(i).seq{1},'K'));
end

% counting positive sample // 181 positive
count_pos=0;
for i=1:sequence_size
    count_pos = count_pos + length(strfind(DB_Pupylation(i).label{1},'1'));
end

% counting negative sample // 2290 negative
count_neg=0;
for i=1:sequence_size
    count_neg = count_neg + length(strfind(DB_Pupylation(i).label{1},'0'));
end

PupData = cell(count_k, 15);
mirror_size = 10;
sequence_loop = 1;
for i=1:sequence_size
	pos_k = strfind(DB_Pupylation(i).seq{1},'K');
	n_terminus = mirror_size + 1;
	c_terminus = length(DB_Pupylation(i).seq{1}) - mirror_size;	

    for j=1:length(pos_k)

        label_sam = str2double(DB_Pupylation(i).label_sep(pos_k(j)));

        PupData{sequence_loop,1} = i;
        PupData{sequence_loop,2} = pos_k(j);
        PupData{sequence_loop,3} = label_sam;

        if (pos_k(j) < n_terminus)
			%ASA
			matrixab = DB_Pupylation(i).ASA(1 : pos_k(j) + mirror_size,:);
			matrixb = DB_Pupylation(i).ASA(pos_k(j)*2 : pos_k(j) + mirror_size,:);
			matrixflip = [flipud(matrixb);matrixab];
			PupData{sequence_loop,4} = matrixflip;
			
			%Phi
			matrixab = DB_Pupylation(i).Phi(1 : pos_k(j) + mirror_size,:);
			matrixb = DB_Pupylation(i).Phi(pos_k(j)*2 : pos_k(j) + mirror_size,:);
			matrixflip = [flipud(matrixb);matrixab];
			PupData{sequence_loop,5} = (matrixflip);
			
			%Psi
			matrixab = DB_Pupylation(i).Psi(1 : pos_k(j) + mirror_size,:);
			matrixb = DB_Pupylation(i).Psi(pos_k(j)*2 : pos_k(j) + mirror_size,:);
			matrixflip = [flipud(matrixb);matrixab];
			PupData{sequence_loop,6} = (matrixflip);
			
			%Theta
			matrixab = DB_Pupylation(i).Theta(1 : pos_k(j) + mirror_size,:);
			matrixb = DB_Pupylation(i).Theta(pos_k(j)*2 : pos_k(j) + mirror_size,:);
			matrixflip = [flipud(matrixb);matrixab];
			PupData{sequence_loop,7} = (matrixflip);
			
			
			%Tau
			matrixab = DB_Pupylation(i).Tau(1 : pos_k(j) + mirror_size,:);
			matrixb = DB_Pupylation(i).Tau(pos_k(j)*2 : pos_k(j) + mirror_size,:);
			matrixflip = [flipud(matrixb);matrixab];
			PupData{sequence_loop,8} = (matrixflip);
			
			%Pc
			matrixab = DB_Pupylation(i).Pc(1 : pos_k(j) + mirror_size,:);
			matrixb = DB_Pupylation(i).Pc(pos_k(j)*2 : pos_k(j) + mirror_size,:);
			matrixflip = [flipud(matrixb);matrixab];
			PupData{sequence_loop,9} = (matrixflip);
			
			%Pe
			matrixab = DB_Pupylation(i).Pe(1 : pos_k(j) + mirror_size,:);
			matrixb = DB_Pupylation(i).Pe(pos_k(j)*2 : pos_k(j) + mirror_size,:);
			matrixflip = [flipud(matrixb);matrixab];
			PupData{sequence_loop,10} = (matrixflip);

			%PSSM (0nly 11 and 12 used)
			matrixab = DB_Pupylation(i).pssm_prob(1 : pos_k(j) + mirror_size,:);
			matrixb = DB_Pupylation(i).pssm_prob(pos_k(j)*2 : pos_k(j) + mirror_size,:);
			matrixflip = [flipud(matrixb);matrixab];
			PupData{sequence_loop,11} = (matrixflip);

			%SSpre
			matrixab = DB_Pupylation(i).SSpre(1 : pos_k(j) + mirror_size,:);
			matrixb = DB_Pupylation(i).SSpre(pos_k(j)*2 : pos_k(j) + mirror_size,:);
			matrixflip = [flipud(matrixb);matrixab];
			PupData{sequence_loop,12} = (matrixflip);

			
        elseif (pos_k(j) >  c_terminus)
			%ASA
			matrixab = DB_Pupylation(i).ASA(pos_k(j) - mirror_size : length(DB_Pupylation(i).seq{1}),:);
			matrixb = DB_Pupylation(i).ASA(pos_k(j) - mirror_size : pos_k(j) - (length(DB_Pupylation(i).seq{1}) - pos_k(j) + 1),:);
			matrixflip = [matrixab; flipud(matrixb)];
			PupData{sequence_loop,4} = matrixflip;

			%Phi
			matrixab = DB_Pupylation(i).Phi(pos_k(j) - mirror_size : length(DB_Pupylation(i).seq{1}),:);
			matrixb = DB_Pupylation(i).Phi(pos_k(j) - mirror_size : pos_k(j) - (length(DB_Pupylation(i).seq{1}) - pos_k(j) + 1),:);
			matrixflip = [matrixab; flipud(matrixb)];
			PupData{sequence_loop,5} = (matrixflip);

			%Psi
			matrixab = DB_Pupylation(i).Psi(pos_k(j) - mirror_size : length(DB_Pupylation(i).seq{1}),:);
			matrixb = DB_Pupylation(i).Psi(pos_k(j) - mirror_size : pos_k(j) - (length(DB_Pupylation(i).seq{1}) - pos_k(j) + 1),:);
			matrixflip = [matrixab; flipud(matrixb)];
			PupData{sequence_loop,6} = (matrixflip);
			
			%Theta
			matrixab = DB_Pupylation(i).Theta(pos_k(j) - mirror_size : length(DB_Pupylation(i).seq{1}),:);
			matrixb = DB_Pupylation(i).Theta(pos_k(j) - mirror_size : pos_k(j) - (length(DB_Pupylation(i).seq{1}) - pos_k(j) + 1),:);
			matrixflip = [matrixab; flipud(matrixb)];
			PupData{sequence_loop,7} = (matrixflip);

			%Tau
			matrixab = DB_Pupylation(i).Tau(pos_k(j) - mirror_size : length(DB_Pupylation(i).seq{1}),:);
			matrixb = DB_Pupylation(i).Tau(pos_k(j) - mirror_size : pos_k(j) - (length(DB_Pupylation(i).seq{1}) - pos_k(j) + 1),:);
			matrixflip = [matrixab; flipud(matrixb)];
			PupData{sequence_loop,8} = (matrixflip);
			
			%Pc
			matrixab = DB_Pupylation(i).Pc(pos_k(j) - mirror_size : length(DB_Pupylation(i).seq{1}),:);
			matrixb = DB_Pupylation(i).Pc(pos_k(j) - mirror_size : pos_k(j) - (length(DB_Pupylation(i).seq{1}) - pos_k(j) + 1),:);
			matrixflip = [matrixab; flipud(matrixb)];
			PupData{sequence_loop,9} = (matrixflip);

			%Pe
			matrixab = DB_Pupylation(i).Pe(pos_k(j) - mirror_size : length(DB_Pupylation(i).seq{1}),:);
			matrixb = DB_Pupylation(i).Pe(pos_k(j) - mirror_size : pos_k(j) - (length(DB_Pupylation(i).seq{1}) - pos_k(j) + 1),:);
			matrixflip = [matrixab; flipud(matrixb)];
			PupData{sequence_loop,10} = (matrixflip);
			
			%PSSM
			matrixab = DB_Pupylation(i).pssm_prob(pos_k(j) - mirror_size : length(DB_Pupylation(i).seq{1}),:);
			matrixb = DB_Pupylation(i).pssm_prob(pos_k(j) - mirror_size : pos_k(j) - (length(DB_Pupylation(i).seq{1}) - pos_k(j) + 1),:);
			matrixflip = [matrixab; flipud(matrixb)];
			PupData{sequence_loop,11} = (matrixflip);		
			
			%SSpre
			matrixab = DB_Pupylation(i).SSpre(pos_k(j) - mirror_size : length(DB_Pupylation(i).seq{1}),:);
			matrixb = DB_Pupylation(i).SSpre(pos_k(j) - mirror_size : pos_k(j) - (length(DB_Pupylation(i).seq{1}) - pos_k(j) + 1),:);
			matrixflip = [matrixab; flipud(matrixb)];
			PupData{sequence_loop,12} = (matrixflip);	
			
		else
			PupData{sequence_loop,4} = (DB_Pupylation(i).ASA(pos_k(j)-mirror_size: pos_k(j)+mirror_size,:)); 
			PupData{sequence_loop,5} = (DB_Pupylation(i).Phi(pos_k(j)-mirror_size: pos_k(j)+mirror_size,:));
			PupData{sequence_loop,6} = (DB_Pupylation(i).Psi(pos_k(j)-mirror_size: pos_k(j)+mirror_size,:));
			PupData{sequence_loop,7} = (DB_Pupylation(i).Theta(pos_k(j)-mirror_size: pos_k(j)+mirror_size,:));
			PupData{sequence_loop,8} = (DB_Pupylation(i).Tau(pos_k(j)-mirror_size: pos_k(j)+mirror_size,:));
			PupData{sequence_loop,9} = (DB_Pupylation(i).Pc(pos_k(j)-mirror_size: pos_k(j)+mirror_size,:));
			PupData{sequence_loop,10} = (DB_Pupylation(i).Pe(pos_k(j)-mirror_size: pos_k(j)+mirror_size,:));
			PupData{sequence_loop,11} = (DB_Pupylation(i).pssm_prob(pos_k(j)-mirror_size: pos_k(j)+mirror_size,:)); %PSSm
			PupData{sequence_loop,12} = (DB_Pupylation(i).SSpre(pos_k(j)-mirror_size: pos_k(j)+mirror_size,:)); %SSPre
		end

		sequence_loop = sequence_loop + 1;
        
    end
end



%--------------------------------------------choose 1------------------------------------------------%
% min - max normalize column 5
%PSSM
newmin = 0;
newmax = 1;
for k=1:count_k

	for m=1:(mirror_size*2)+1
	
		minval = min(PupData{k,11}(m,:)); % to be caluclated in the column
		maxval = max(PupData{k,11}(m,:)); % to be caluclated in the column
			
		for v=1:length(PupData{k,11}(m,:))
		
			val = PupData{k,11}(m,v); % to be inserted
			valnew = (val - minval)/(maxval-minval)* (newmax - newmin) + newmin;
			PupData{k,16}(m,v) = valnew;
		
		end
	end
end


%--------------------------------------------choose 1------------------------------------------------%
%zscore normalization
% for k=1:count_k


	% for m=1:(mirror_size*2)+1

		% avg_val = mean(PupData{k,11}(m,:));
		% st_val = std(PupData{k,11}(m,:));
	
		% for v=1:length(PupData{k,11}(m,:))
		
			% val = PupData{k,11}(m,v); % to be inserted
			% valnew = val - (avg_val/st_val);
			% PupData{k,17}(m,v) = valnew;
		
		% end
	% end
% end

%--------------------------------------------choose 1------------------------------------------------%

%bigram on PSSM

for z=1:count_k
    PSSM_initial_matrix = PupData{z,16};
    sizepssm_i = size(PSSM_initial_matrix);
    bigram_L = sizepssm_i(1);
    bigram_row_col = sizepssm_i(2);
    bigram_matrix = [];

    for bigram_row=1:bigram_row_col
        for bigram_col=1:bigram_row_col
            bigram_ij_sum_vector = 0;

            for i=1:bigram_L - 1		
                bigram_ij_sum_vector = bigram_ij_sum_vector + (PSSM_initial_matrix(i,bigram_row) * PSSM_initial_matrix(i+1,bigram_col));
            end

            bigram_matrix(bigram_row, bigram_col) = bigram_ij_sum_vector;
        end
    end

    bigram_matrix_transpose = bigram_matrix'; 
    bigram_vector = bigram_matrix_transpose(:)';
	PupData{z,11} = bigram_vector;
end

% %bigram on SSPre
% for z=1:count_k
    % PSSM_initial_matrix = PupData{z,12};
    % sizepssm_i = size(PSSM_initial_matrix);
    % bigram_L = sizepssm_i(1);
    % bigram_row_col = sizepssm_i(2);
    % bigram_matrix = [];

    % for bigram_row=1:bigram_row_col
        % for bigram_col=1:bigram_row_col
            % bigram_ij_sum_vector = 0;

            % for i=1:bigram_L - 1		
                % bigram_ij_sum_vector = bigram_ij_sum_vector + (PSSM_initial_matrix(i,bigram_row) * PSSM_initial_matrix(i+1,bigram_col));
            % end

            % bigram_matrix(bigram_row, bigram_col) = bigram_ij_sum_vector;
        % end
    % end

    % bigram_matrix_transpose = bigram_matrix'; 
    % bigram_vector = bigram_matrix_transpose(:)';
	% PupData{z,12} = bigram_vector;
% end

%PSSM & SSprr
for k=1:count_k
   % PupData{k,13} = [PupData{k,11} PupData{k,12} ];
    PupData{k,13} = [PupData{k,11} ]; % just PSSM
	featurevector = PupData{k,13}';
	PupData{k,14} = featurevector(:)';
	%PupData{k,14} = [featurevector(:)' PupData{k,12}];
end

% selected or not
summax = 23;

for z=1:count_k
	DistData_single =[];
	
	if(PupData{z,3} == 1)
		PupData{z,15} = 1;
		PupData{z,16} = 1;
	else
	
		for j=1:count_k
					
			% dist = sum(sum(pdist2(PupData{z,4}, PupData{j,4},'minkowski')));
		
            dist = pdist2(PupData{z,14}, PupData{j,14},'minkowski');
			DistData_single = [DistData_single; dist, PupData{j,3}];
		end
		
		sorted_DistData = sortrows(DistData_single, 1);
		sumzero = sum(sorted_DistData(1:summax,2));
		
		PupData{z,16} = sorted_DistData;
		
		if (sumzero >= 1)
			PupData{z,15} = 0;
		else
			PupData{z,15} = 1;
		end
	end
end
%kNN
%k
%22 = 361 samples (181 positive and 180 negative) ratio (positive:negative) -> 1:1
%21 = 400 samples (181 positive and 219 negative) ratio (positive:negative) -> 1:1.21
%20 = 448 samples (181 positive and 267 negative) ratio (positive:negative) -> 1:1.48
%19 = 494 samples (181 positive and 313 negative) ratio (positive:negative) -> 1:1.73
%18 = 572 samples (181 positive and 391 negative) ratio (positive:negative) -> 1:2.16
% 30 max

%=====run from here===============
count_k=2471;
summax = 70;

for j=1:length(PupData)
    
    if(PupData{j,3} == 1)
		PupData{j,15} = 1;
    else
        
		sumzero = sum(PupData{j,16}(1:summax,2));
		
		if (sumzero >= 1)
			PupData{j,15} = 0;
		else
			PupData{j,15} = 1;
        end
    end
end
	
	
s=0;
for j=1:length(PupData)
s  = s + PupData{j,15};

end
s

%sample data
postive_s = 181;
negative_s = s - 181;

postive_s_i = 1;
negative_s_i = 1;
s_i = 1;

SamplePupData = cell(s, 5);
PostiveSamplePupData = cell(postive_s, 5);
NegativeSamplePupData = cell(negative_s, 5);

% for j=1:count_k
    
    % if(PupData{j,15} == 1)
		% SamplePupData(s_i,1:3) = PupData(j,1:3);
        % %for compare with other result, index is needed
		% SamplePupData{s_i,4} = j;
		% SamplePupData{s_i,5} = PupData{j,14};
		
        % if (PupData{j,3} == 1)
			% PostiveSamplePupData(postive_s_i,1:3) = PupData(j,1:3);
			% %for compare with other result, index is needed
			% PostiveSamplePupData{postive_s_i,4} = j;
			% PostiveSamplePupData{postive_s_i,5} = PupData{j,14};
            % postive_s_i = postive_s_i + 1;
		% else
			% NegativeSamplePupData(negative_s_i,1:3) = PupData(j,1:3);
			% %for compare with other result, index is needed
			% NegativeSamplePupData{negative_s_i,4} = j;
			% NegativeSamplePupData{negative_s_i,5} = PupData{j,14};
            % negative_s_i = negative_s_i + 1;
        % end
        
        % s_i = s_i + 1;
    % end
% end


for j=1:count_k
    
    if(PupData{j,15} == 1)
		SamplePupData(s_i,1:3) = PupData(j,1:3);
        %for compare with other result, index is needed
		SamplePupData{s_i,4} = j;
		SamplePupData{s_i,5} = PupData{j,14};
		
        if (PupData{j,3} == 1)
			PostiveSamplePupData(postive_s_i,1:3) = PupData(j,1:3);
			%for compare with other result, index is needed
			PostiveSamplePupData{postive_s_i,4} = j;
			PostiveSamplePupData{postive_s_i,5} = PupData{j,14};
            postive_s_i = postive_s_i + 1;
		else
			NegativeSamplePupData(negative_s_i,1:3) = PupData(j,1:3);
			%for compare with other result, index is needed
			NegativeSamplePupData{negative_s_i,4} = j;
			NegativeSamplePupData{negative_s_i,5} = PupData{j,14};
            negative_s_i = negative_s_i + 1;
        end
        
        s_i = s_i + 1;
    end
end


% k fold
foldsvalue = [6 8 10];
foldResult = [];
ComparefoldResult = [];
resultindex=1; % fill in the result matlab

for foldsloop = 1:3
	k = foldsvalue(foldsloop);
	k_p = PostiveSamplePupData;
	k_n = NegativeSamplePupData;

	pos_kfold = kfolddiv( k, k_p );
	neg_kfold = kfolddiv( k, k_n );


	TestData = cell(k, 1);

	for i=1:k
		foldpositive = PostiveSamplePupData(pos_kfold{i,3}(1):pos_kfold{i,3}(2), 5);
		foldpositive(:,2) = PostiveSamplePupData(pos_kfold{i,3}(1):pos_kfold{i,3}(2), 3);
		%for compare with other result, index is needed
		foldpositive(:,3) = PostiveSamplePupData(pos_kfold{i,3}(1):pos_kfold{i,3}(2), 4);

		foldnegative = NegativeSamplePupData(neg_kfold{i,3}(1):neg_kfold{i,3}(2), 5);
		foldnegative(:,2) = NegativeSamplePupData(neg_kfold{i,3}(1):neg_kfold{i,3}(2), 3);
		%for compare with other result, index is needed
		foldnegative(:,3) = NegativeSamplePupData(neg_kfold{i,3}(1):neg_kfold{i,3}(2), 4);

		testfold = [foldpositive; foldnegative];
		TestData{i,1} = testfold;

	end

	if (k == 6)
		TrainData = cell(k, 1);
		TrainData{1,1} = [TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}];
		TrainData{2,1} = [TestData{1,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}];
		TrainData{3,1} = [TestData{1,1}; TestData{2,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}];
		TrainData{4,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{5,1}; TestData{6,1}];
		TrainData{5,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{6,1}];
		TrainData{6,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}];
		

	end

	%8 folds
	if (k == 8)
		TrainData = cell(k, 1);
		TrainData{1,1} = [TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}; TestData{7,1}; TestData{8,1}];
		TrainData{2,1} = [TestData{1,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}; TestData{7,1}; TestData{8,1}];
		TrainData{3,1} = [TestData{1,1}; TestData{2,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}; TestData{7,1}; TestData{8,1}];
		TrainData{4,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{5,1}; TestData{6,1}; TestData{7,1}; TestData{8,1}];
		TrainData{5,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{6,1}; TestData{7,1}; TestData{8,1}];
		TrainData{6,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{7,1}; TestData{8,1}];
		TrainData{7,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}; TestData{8,1}];
		TrainData{8,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}; TestData{7,1}];


	end

	%10 folds
	if (k == 10)
		TrainData = cell(k, 1);
		TrainData{1,1} = [TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}; TestData{7,1}; TestData{8,1}; TestData{9,1}; TestData{10,1}];
		TrainData{2,1} = [TestData{1,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}; TestData{7,1}; TestData{8,1}; TestData{9,1}; TestData{10,1}];
		TrainData{3,1} = [TestData{1,1}; TestData{2,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}; TestData{7,1}; TestData{8,1}; TestData{9,1}; TestData{10,1}];
		TrainData{4,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{5,1}; TestData{6,1}; TestData{7,1}; TestData{8,1}; TestData{9,1}; TestData{10,1}];
		TrainData{5,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{6,1}; TestData{7,1}; TestData{8,1}; TestData{9,1}; TestData{10,1}];
		TrainData{6,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{7,1}; TestData{8,1}; TestData{9,1}; TestData{10,1}];
		TrainData{7,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}; TestData{8,1}; TestData{9,1}; TestData{10,1}];
		TrainData{8,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}; TestData{7,1}; TestData{9,1}; TestData{10,1}];
		TrainData{9,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}; TestData{7,1}; TestData{8,1}; TestData{10,1}];
		TrainData{10,1} = [TestData{1,1}; TestData{2,1}; TestData{3,1}; TestData{4,1}; TestData{5,1}; TestData{6,1}; TestData{7,1}; TestData{8,1}; TestData{9,1}];


	
	end

	initial_index = resultindex;
	%train and test using Matlab
	for i=1:k
		traindata = TrainData{i,1};
		testdata = TestData{i,1};
		
		foldResult(resultindex,:) = PerformanceAssessment( traindata, testdata);
		
		% comapre
		sizetestdata = length(testdata);
		comparetestdata = [];
		for compare_i=1:sizetestdata
			if (cc{testdata{compare_i,3},3} == 1)
				comparetestdata(compare_i,1) = 1;
			else
				comparetestdata(compare_i,1) = 0;
			end
		
		end
		
		ComparefoldResult(resultindex,:) = ComparePerformanceAssessment( testdata, comparetestdata);
		
		resultindex = resultindex + 1;
	end

	finalindex = resultindex-1;
	for j=5:9
		foldResult(resultindex, j) = mean(foldResult(initial_index:finalindex, j));
		
		% compare
		ComparefoldResult(resultindex, j) = mean(ComparefoldResult(initial_index:finalindex, j));
	end

	resultindex = resultindex + 2;

end