N = 3005;

cell_class_N = cell(N,1);
cell_group_num_N = cell(N,1);
cell_tissue_N = cell(N,1);
cell_subclass_N = cell(N,1);

for i = 1:N
    ind = find(strcmp(cell_id, cell_id_ordered{i}));
    cell_class_N{i} = cell_class{ind};
    cell_group_num_N{i} = cell_group_num{ind};
    cell_tissue_N{i} = cell_tissue{ind};
    cell_subclass_N{i} = cell_subclass{ind};
end

%% 

save_to_label_file('cell_class_3005.txt',cell_class_N,'./metadata');
save_to_label_file('cell_group_num_3005.txt',cell_group_num_N,'./metadata');
save_to_label_file('cell_tissue_3005.txt',cell_tissue_N,'./metadata');
save_to_label_file('cell_subclass_3005.txt',cell_subclass_N,'./metadata');