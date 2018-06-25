%% Emulated 4-PAM decoding

function [ pbits ] = dec_4pam( psyms )

    L = numel(psyms);
    pbits = zeros(2*L,1);
    for l = 1:L
        switch psyms(l)
            case -3
                pbits(2*l-1:2*l) = [0 ; 0];
            case -1
                pbits(2*l-1:2*l) = [0 ; 1];
            case 1
                pbits(2*l-1:2*l) = [1 ; 1];
            case 3
                pbits(2*l-1:2*l) = [1 ; 0];
        end
    end

end