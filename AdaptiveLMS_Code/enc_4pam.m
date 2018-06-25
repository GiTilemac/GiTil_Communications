%% Emulated 4-PAM Encoding

function [ psyms ] = enc_4pam( pbits )

    L = numel(pbits);
    psyms = zeros(L/2,1);
    for l = 1:2:L-1
        enc_pair = strcat(num2str(pbits(l)),num2str(pbits(l+1)));
        switch enc_pair
            case '00'
                psyms(ceil(l/2)) = -3;
            case '01'
                psyms(ceil(l/2)) = -1;
            case '11'
                psyms(ceil(l/2)) = 1;
            case '10'
                psyms(ceil(l/2)) = 3;
        end   
        
    end
    
end