%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  subpanel.m: Identifying substitutable gene sets 
%%  and generating all possible combinations of them 
%%  as potential biomarker panels for a given GOBP result
%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. import 82 genes
fid = fopen('Gene82.txt', 'r');
data = textscan(fid, '%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'headerlines', 1);

gene82 = data{1,1};
gene82alias = data{1,2};
gene11ids = data{1,3};
gene11 = gene82(gene11ids == 1);
gene11alias = gene82alias(gene11ids == 1);
gene82fold = data{1,5};
gene11fold = gene82fold(gene11ids == 1);
clear fid data gene11ids;

%% 2. import GOBPs
data = readtable('GOBP_82genes_alias.txt', ...
    'delimiter', '\t');

GOBP82.id = cell(length(data.Term), 1);
GOBP82.term = cell(length(data.Term), 1);
GOBP82.gene = struct;
for i = 1:length(data.Term)
    tm = strsplit(data.Term{i}, '~');
    GOBP82.id(i, 1) = tm(1,1);
    GOBP82.term(i, 1) = tm(1, 2);
    
    tm = strsplit(data.Genes{i}, ', ');
    GOBP82.gene(i, 1).symbol = tm';
end

GOBP82.count = data.Count;
GOBP82.pVal = data.PValue;
clear tm data;


%% 3. Selection of GOBPs including 11 genes
selids = [];
for i = 1:length(GOBP82.gene)
    tid = find(ismember(GOBP82.gene(i).symbol, gene11alias) == 1);
    if ~isempty(tid) == 1
        selids = [selids; i];
    end
end

selids = setdiff([1:1:503], selids)';


%% 4. List of genes for substitution
geneSubAlias = struct;
geneSub = struct;
for i = 1:length(gene11alias)
    geneSubAlias(i, 1).list = [];
    
    for j = 1:length(selids)
        tGenes = GOBP82.gene(selids(j)).symbol;
        tid = find(ismember(tGenes, gene11alias(i)) == 1);
        if ~isempty(tid) == 1 && length(tGenes) > 1
            geneSubAlias(i, 1).list = [geneSubAlias(i, 1).list; tGenes];
        end
    end
    
    if ~isempty(geneSubAlias(i, 1).list) == 1
        geneSubAlias(i, 1).list = unique(geneSubAlias(i, 1).list);
        geneSubAlias(i, 1).list = setdiff(geneSubAlias(i, 1).list, gene11alias);
        
        if gene11fold(i) > 0
            retainGenes = gene82alias(gene82fold > 0);
            geneSubAlias(i, 1).list = intersect(geneSubAlias(i, 1).list, retainGenes);
        elseif gene11fold(i) < 0
            retainGenes = gene82alias(gene82fold < 0);
            geneSubAlias(i, 1).list = intersect(geneSubAlias(i, 1).list, retainGenes);
        end
        
        geneSub(i, 1).list = gene82(ismember(gene82alias, geneSubAlias(i, 1).list) == 1);
    end
    
end