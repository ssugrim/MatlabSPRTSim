[chan,lam,del,split] = size(Data_SPRT_channels_found_H0);

err_2 = 0.05;
err_1 = 1;

max_disc = 0;

for c = 1:chan
    for l = 1:lam
        for d = 1:del
            for s = 1:split
                if Data_SPRT_error_1(c,l,d,s) < err_1 && Data_SPRT_error_2(c,l,d,s) < err_2
                    if Data_SPRT_channels_found_H0(c,l,d,s) > max_disc
                        max_disc = Data_SPRT_channels_found_H0(c,l,d,s);
                        sprintf('Max:%5.2f at Channels:%5.2f, Lambda:%5.2f, delta:%5.2f, split:%5.2f, error_1:%5.2f, error_2:%5.2f',max_disc,Channels(c),lambda_range(l),delta_range(d),split_range(s),Data_SPRT_error_1(c,l,d,s),Data_SPRT_error_2(c,l,d,s))
                    end                    
                end
            end
        end
    end
end


