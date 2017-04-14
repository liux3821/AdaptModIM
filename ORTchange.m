function [ort_change,delta_ort,old_ort] = ORTchange(taskRandomOrt,dpp,xCenter,w,ort_change,delta_ort,old_ort)

[adj_x,adj_y,buttons] = GetMouse(w);
if adj_x > xCenter
    if delta_ort > 0
        if ort_change > old_ort
            old_ort = floor(ort_change);
        else
            old_ort = ceil(ort_change);
        end
    end
    ort_change = ((adj_x - xCenter)*dpp)+taskRandomOrt;
    delta_ort = floor(abs(old_ort-ort_change));
else
    if delta_ort > 0
        if ort_change > old_ort
            old_ort = floor(ort_change);
        else
            old_ort = ceil(ort_change);
        end
    end
    ort_change = taskRandomOrt-((xCenter - adj_x)*dpp);
    delta_ort = floor(abs(ort_change - old_ort));
end

