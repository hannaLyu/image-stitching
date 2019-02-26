function varargout = brief_matching(descriptor1, descriptor2)
    matching_mode = 'burteforce';% kd tree version will be provided
    if strcmp(matching_mode,'burteforce')
        n1 = size(descriptor1,1);
        n2 = size(descriptor2,1);
        matching_pairs = zeros(n1,2);
        for i = 1:n1
            descriptor1rep = repmat(descriptor1(i,:),n2,1);
            despxor = xor(descriptor1rep,descriptor2);
            despxor = sum(despxor,2);
            [~,id] = min(despxor);
%             [val,minid] = sort(despxor,'ascend');
%             best = val(1);
%             secondbest = val(2);
%             if best < secondbest*0.9
                matching_pairs(i,:) = [i, id];
%             else
%                 matching_pairs(i,:) = [i, inf];
%             end
        end
    end
    varargout{1} = matching_pairs;
end