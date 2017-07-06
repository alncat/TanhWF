count = 0;
for i = 1:19
    Example_1D_Gaussian;
    if min(Relerrs) < 0.01
        count = count + 1;
    end
end
display(count);