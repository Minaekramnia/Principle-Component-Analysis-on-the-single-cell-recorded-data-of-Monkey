%Main Goal : PCA on the data

%First Step : preparing the main matrix that we want apply PCA . which is
%firing rate for each bin in each condition (ce,Pd,Respdir) . its dim will
%be nc*t


%list folders of cells in an array(cause they are not in order for +1 adding e.x 0095)
%go inside each folders with a for loop
%load the arbitrary .mat file that we want to use. using sprintf [C1 not
%C0???] ** I should know about the i cell  from last step because I need
% its name for finding the mat file.
%check if exist load that if else continue and go to other cells
%write some functions for finding Trials with the same condition :


cd I:\lists
mylist= textread ('V1list.txt', '%s');
s1 = '/bgc';
s2 = 'I:';
mylist= strrep(mylist, s1, s2);
mylist= strrep(mylist, '/','\');
nFile = numel(mylist);

T= 4; % minimum number of Trials.

for fileNum = 1 : nFile
    a=strfind(mylist{fileNum}(2:end-1),'PCEDP');
    p = [];
    c = [];
    ch = [];
    if ~isempty(a)
        load (mylist{fileNum}(2:end-1))
        z = Expt.Trials;
        if (length(z) > T)
            for i = 1 : length (z)
                if (~isempty(z(i).Pd) && ~isempty(z(i).ce) && ~isempty(z(i).RespDir))
                    p(i) = z(i).Pd; % matrix of Pd for each cell
                    c(i) = z(i).ce; % matrix of ce for each cell
                    ch(i) = z(i).RespDir; % matrix of choices for each cell
                end
            end
            pdu= unique (p);
            if size(pdu) == [1,2]
                cu = unique (c);
                chu = unique (ch);
                count = 1;
                for l = 1 : length(cu)
                    for j= 1 : length(pdu)
                        for k = 1 : length(chu)
                            if(chu(k)~=0)
                                if (0 < cu(l) && cu(l) <= 0.25)
                                    ind = find ( round (10000*c)== round (10000 *cu(l)) & round (10000 * p) == round (10000* pdu(j))  & ch == chu(k)); % index of trials with the same condition
                                    if (size(ind, 2) > T)
                                        cond{fileNum}(count, :) = [cu(l) pdu(j) chu(k)];

                                        for m = 1 : length (ind)
                                            spikesPerBin(m, :) = histc(z(ind(m)).Spikes,0 : 1000 : 20000) ; %100msec
                                        end
                                        ave = mean(spikesPerBin, 1);%average spikes of all trials in the same cond
                                        fRate{fileNum, 1}(count, :) = ave * 10 ; % Firing Rate(Sp/sec) in each bin per condition
                                        count = count + 1;
                                    end
                                elseif(cu(l)==0)
%                                     disp('heeeeeeeeeeey')
                                    ind = find ( round (10000*c)== 0 & ch == chu(k)); % index of trials with the same condition
                                    if (size(ind, 2) > T)
                                        cond{fileNum}(count, :) = [cu(l) 0 chu(k)]; %is it right?

                                        for m = 1 : length (ind)
                                            spikesPerBin(m, :) = histc(z(ind(m)).Spikes,0 : 1000 : 20000) ;
                                        end
                                        ave = mean(spikesPerBin, 1);
                                        fRate{fileNum, 1}(count, :) = ave * 10;  % Firing Rate(S/sec) in each bin per condition
                                        count = count + 1;


                                    end
                                end
                            end
                        end


                    end


                end
                % disp(fRate{fileNum}(:, :))
                
            end
        end

    end

end
%fRateMat = cell2mat(fRate);
%[COEFF,SCORE]=pca(fRateMat')