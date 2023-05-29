% clear all, clc;
%%  
load('PNP_name.mat');  load('stream_85.mat', 'down', 'up');
c = containers.Map('KeyType','char','ValueType','double');
% Protein and gene name are identical
for i = 1:length(down)
    if c.isKey(PNP_down_name(i))
        c(PNP_down_name(i)) = (c(PNP_down_name(i)) + down(i)) / 2;
    else
        c(PNP_down_name(i)) = down(i);
    end
end

for i = 1:length(up)
    if c.isKey(PNP_up_name(i))
        c(PNP_up_name(i)) = (c(PNP_up_name(i)) + up(i)) / 2;
    else
        c(PNP_up_name(i)) = up(i);
    end
end


keys = c.keys;
values = cell2mat(c.values);
[rank_value, sortId] = sort(values,'descend');
rank_value = rank_value';
rank_gene = cell(length(sortId),1);
for i = 1:length(sortId)
    rank_gene{i} = keys{sortId(i)};
end
rank_gene = string(rank_gene);
save('Rank.mat','rank_value', 'rank_gene');

%% 
File = fopen('rank.txt','a');
for i = 1:length(rank_gene)
    if i == length(rank_gene)
        fprintf(File,'%s',rank_gene(i));
    else
        fprintf(File,'%s\n',rank_gene(i));
    end
end
fclose(File);
