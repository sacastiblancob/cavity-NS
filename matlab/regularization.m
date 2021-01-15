function [DUpr] = regularization(DUp,MLap)
    [U] = svds(MLap,1,'smallest');
    DUpr = (1-U.^2).*DUp;
end