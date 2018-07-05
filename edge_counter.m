clear all
close all
clc

load dataset2.mat

k_end = 2000;

edge_counts = zeros(17,1);
for i=1:length(l)
    edge_counts(i) = nnz(r(1:k_end,i));
end