%%
FASTAfilename = 'D:\human_g1k_v37.fasta';
fileInfo = dir(which(FASTAfilename));
fidIn = fopen(FASTAfilename,'r');
header = fgetl(fidIn);
fidOut = fopen('human_g1k_v37_prepared.fasta','at');

% Removing the spaces and make the whole Genome in one line
newLine = sprintf('\n');
blockSize = 2^21;
tic
while ~feof(fidIn)
    % Read in the data
    GenomeData = fread(fidIn,blockSize,'*char')';
    % Remove new lines
    GenomeData = strrep(GenomeData,newLine,'');
    % Save the Whole Genome in New File
    fprintf(fidOut, '%s', GenomeData);
end
toc
fclose(fidIn);
fclose(fidOut);

%%
fidIn = fopen('human_g1k_v37_prepared.fasta','r');
fidOut = fopen('out.fastq', 'w');
tic
for i = 1:1000000000
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % mrFAST Alignment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Structure = samread('D:\ERR240727_1.map', 'blockread', i);
    Read = Structure.Sequence;
    % Read in the Reference
    status = fseek(fidIn,Structure.Position-1,'bof');
    Ref = fread(fidIn, 100,'*char')';
    Cigar = Structure.CigarString;

    fprintf(fidOut, '%s %s %s\n', Read, Ref, Cigar);
end
toc
fclose(fidIn);
fclose(fidOut);