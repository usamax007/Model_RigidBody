function f_out = rescaleSinusoid(f_in,f_mean,f_amp)
% amp is peak to peak    
f_temp = f_in - mean(f_in);
if((max(f_temp) - min(f_temp)) == 0)
    f_out = f_mean*ones(size(f_in));
else
    f_temp = f_temp/(max(f_temp) - min(f_temp));
    f_out = f_mean + f_temp*f_amp;
end

end