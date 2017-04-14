function [old_ort,ort_change,delta_ort] = ort_change(ort_change,adj_x,xCenter,dpp,taskRandomOrt)

old_ort = ort_change;
ort_change = ((adj_x - xCenter)*dpp)+taskRandomOrt;
delta_ort = round(abs(ort_change - old_ort));