function neighbors = find_neighbors_jhs(target,tri,radius)
for i=1:radius
    if i ==1
        neighbors = unique(spm_vec(tri(mod(find(tri == target),size(tri,1)),:)));
    else
        neighbors_new = [];
        for j = 1:length(neighbors)
            neighbors_new = union(neighbors_new,unique(spm_vec(tri(mod(find(tri == neighbors(j)),size(tri,1)),:))));
        end
        neighbors = neighbors_new;
    end
end
end