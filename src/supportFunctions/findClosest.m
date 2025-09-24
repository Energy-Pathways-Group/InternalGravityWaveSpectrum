% Function to find the closest value in kRadial
function idx = findClosest(kRadial, k)
    [~, idx] = min(abs(kRadial - k));  % Finds the index of the closest value in kRadial to k
end