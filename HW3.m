%HW3

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'

% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 
matchval=2;
mismatchval=-1;
ofdiag=ones(4)-eye(4);
S=matchval*eye(4)+mismatchval*ofdiag;
seq1='GTAATCC';
seq2='GTATCCG';
figure;
[score,align,start]=swalign(seq1,seq2,'Alphabet','nt','ScoringMatrix',S,'Gapopen',1,'Showscore',true);

%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 
accession='NM_002746';
gb_data=getgenbank(accession);
seq1=gb_data.Sequence;

accession2='NM_002745';
gb_data2=getgenbank(accession2);
seq2=gb_data2.Sequence;

[score,align,start]=swalign(seq1,seq2,'Alphabet','nt','Showscore',true);


%Using ncbi blastn tools result: 790bps, 41%.

% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?
NP_002737

NP_002736

accession='NM_002746';
gb_data=getgenbank(accession);
seq1=gb_data.Sequence;
aa1=gb_data.CDS.translation;

accession2='NM_002745';
gb_data2=getgenbank(accession2);
aa2=gb_data.CDS.translation;

[score,align,start]=swalign(aa1,aa2,'Alphabet','nt','Showscore',true);

%Using ncbi blastp tools result: 305/346; 88%

% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 
X64605
D10939
accession='X64605';
gb_data=getgenbank(accession);
seq1=gb_data.Sequence;

accession2='D10939';
gb_data2=getgenbank(accession2);
seq2=gb_data2.Sequence;

[score,align,start]=swalign(aa1,aa2,'Alphabet','nt','Showscore',true);
% Using ncbi blastn tools result: 77% identity.
%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 
function xx=blastN(a,b)
aa=getgenbank(a);
seq=aa.Sequence;
[requestID,requestTime]=blastncbi(seq,'blastn');
blast_data=getblast2(requestID,'WaitTime',requestTime);
ax=[];
for i=1:b;
    ax(i)=blastdata.Hits(i).Name;
end
xx=ax
end

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 
function F=blastNmatch(a)
aaa={};
bb={};
aa=getgenbank(a);
seq=aa.Sequence;
[requestID,requestTime]=blastncbi(seq,'blastn');
blast_data=getblast(requestID,'WaitTime',requestTime);
N=length(blast_data.Hits);
cc=[];
dd=[];
for i=1:N
    bb{i}={blast_data.Hits(i).Name};
    if ~contains(char(bb{i}),'Homo sapiens')
      cc=[cc,i];
    else 
     dd=[dd,i];
    end
end
if isempty(dd)
    aaa=bb{cc(1)};
    fprintf('no human is found',char(aaa{1}));
else if isempty(cc)
     aaa=bb{cc(1)};
    fprintf('no nonhuman is found',char(aaa{1}));
    else
    aaa={char(bb{cc(1)});char(bb{dd(1)})};
    fprintf('the closest human is %s, the closest nonhuman is %s',char(aaa{2}),char(aaa{1}));
    end
end
end
% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 
blastNmatch('NM_002746')
Blast results are not available yet. Please wait ...
the closest human is gi|91718898|ref|NM_002746.2| Homo sapiens mitogen-activated protein kinase 3 (MAPK3), transcript variant 1, mRNA, the closest nonhuman is gi|31220|emb|X60188.1| Human ERK1 mRNA for protein serine/threonine kinase>> 

blastNmatch('NM_002745')
Blast results are not available yet. Please wait ...
the closest human is gi|75709178|ref|NM_002745.4| Homo sapiens mitogen-activated protein kinase 1 (MAPK1), transcript variant 1, mRNA, the closest nonhuman is gi|1034143019|ref|XM_003317123.4| PREDICTED: Pan troglodytes mitogen-activated protein kinase 1 (MAPK1), mRNA>>

%My code focuses on characterizing Homosapiens and other species. However, in the first result, my code seems to fail. I think it is because the criteria of human gene does not limit to 'Homosapiens'.  
