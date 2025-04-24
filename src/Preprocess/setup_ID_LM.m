function model = setup_ID_LM(model)

    for ii = 1:model.nnp
        for jj = 1:model.ndof
            qq = model.ndof*(ii - 1) + jj;
            model.ID(qq) = qq;
            if model.flags(qq) == 2
                model.d(qq) = model.e_bc(qq);
            end
        end
    end

    for ee = 1:model.nel
        for aa = 1:model.nen
            for ii = 1:model.ndof
                pp = model.ndof * (aa - 1) + ii;
                model.LM(pp, ee) = model.ID(ii, model.IEN(aa, ee));
            end
        end
    end


end