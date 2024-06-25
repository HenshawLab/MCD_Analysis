%%
clearvars -except OutputMain ExperimentOutDir ExpName NPos BioReps Reps
mkdir([OutputMain 'Spreadsheets/']);

load([ExperimentOutDir ExpName '_beta_analysis.mat']);
delete([OutputMain 'Spreadsheets/' ExpName '.xls']);

%%

for ip = 1:NPos

    t = beta_output(ip).time';
    beta = beta_output(ip).beta';

    NCols = 2*size(t,2);
    jt = 1:2:NCols; jb = 2:2:NCols;

    datamat = zeros(size(beta,1),NCols);

    datamat(:,jt) = t; 
    datamat(:,jb) = beta;

    bionums = beta_output(ip).BioRep; technums = beta_output(ip).TechRep;

    for i = 1:(size(t,2))
        colname{2*i-1} = ['time_bio-' num2str(BioReps(bionums(i))) '_rep-' num2str(Reps(technums(i)))];
        colname{2*i} = ['beta_bio-' num2str(BioReps(bionums(i))) '_rep-' num2str(Reps(technums(i)))];
    end

    T = array2table(datamat,'VariableNames',colname);
    writetable(T,[OutputMain 'Spreadsheets/' ExpName '.xls'],...
        'Sheet',['Channel' num2str(ip)]);

end