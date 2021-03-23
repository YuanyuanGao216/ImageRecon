function [brain_vertices_new] = down_sample_vertices(brain_vertices, threshod)
brain_vertices_new = [];

for i = 1:length(brain_vertices)
    point = brain_vertices(i,:);
    
    if isempty(brain_vertices_new) % add the first point to the new set
        brain_vertices_new(end+1,:) = point;
        continue
    else
        % is it close to any of the points in new_vertices?
        close_flag = 0;
        for j = 1:size(brain_vertices_new,1)
            distance = sqrt((point(1) - brain_vertices_new(j,1))^2 + ...
                (point(2) - brain_vertices_new(j,2))^2 + ...
                (point(3) - brain_vertices_new(j,3))^2);
            if distance < threshod
                close_flag = 1;
                break
            end
        end
    end
    % after running all the points in the new vertices, there is no close
    % points, we can add this one
    if close_flag == 0
        brain_vertices_new(end + 1,:) = point; 
    end
end