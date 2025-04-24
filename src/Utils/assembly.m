function model = assembly(model, k_ele, f_ele, ee)

    % global assembly
    for aa = 1 : model.nen
        for ii = 1 : model.ndof
            pp = model.ndof * (aa - 1) + ii;
            PP = model.LM(pp, ee);
            model.f(PP) = model.f(PP) + f_ele(pp);
            for bb = 1 : model.nen
                for jj = 1 : model.ndof
                    qq = model.ndof * (bb - 1) + jj;
                    QQ = model.LM(qq, ee);
                    model.K(PP, QQ) = model.K(PP, QQ) + k_ele(pp, qq);
                end % end of jj-loop
            end   % end of for bb
        end % end of ii-loop
    end    % end of for aa

end