% This function computes the modified matrix exponential
% for a given parameter beta.

function [exp_0] = exmpp_0(Matrix, beta)
    exp_0 = expm(beta*Matrix) - eye(size(Matrix));
end
