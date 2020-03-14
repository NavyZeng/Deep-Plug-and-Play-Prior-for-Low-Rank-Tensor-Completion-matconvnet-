function [X] = myNormalization(label)
sz=size(label);

if max(label)>1
    for s = 1:sz(3)
        maxXs = max(max(label(:,:,s)));
        X(:,:,s) =  label(:,:,s)/maxXs;
    end
else
    X = label;
end
end

