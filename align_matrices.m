function [aligned_data1, aligned_data2] = align_matrices(data1, data2)
    m = length(data1);
    n = length(data2);

    if m > n
        delete_count = m - n;
        aligned_data1 = data1(delete_count+1:end);
        aligned_data2 = data2;
    elseif n > m
        delete_count = n - m;
        aligned_data1 = data1;
        aligned_data2 = data2(delete_count+1:end);
    else
        aligned_data1 = data1;
        aligned_data2 = data2;
    end
end