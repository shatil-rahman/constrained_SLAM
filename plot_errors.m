errors_var = [];

for i=1:length(test(:,1))
    if test(i,4) < 3
        errors_var = [errors_var; test(i,4)];
    end
end