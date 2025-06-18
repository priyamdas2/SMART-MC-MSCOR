% Converts treatment sequence matrix with 2 columns into a cell with
% K (no. of patients) indivisual treatment sequences

function output = Trt_Seq_Mat2Cell(TrtSeqMat, N)
SeqTemp = cell(1, N);
PatientNum = 1;
SeqTemp{1} = TrtSeqMat(1, 2);

for i = 2:size(TrtSeqMat, 1)
    if(TrtSeqMat(i,1) == TrtSeqMat(i-1,1))
        SeqTemp{PatientNum} = [SeqTemp{PatientNum} TrtSeqMat(i,2)];
    else
        PatientNum = PatientNum + 1;
        SeqTemp{PatientNum} = TrtSeqMat(i,2);
    end
end

output = SeqTemp;

