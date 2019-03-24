function varargout = bruteforce(descriptor1, descriptor2)

        n1 = size(descriptor1,1);
        n2 = size(descriptor2,1);
        matching_pairs = zeros(n1,2);
        for i = 1:n1
            descriptor1rep = repmat(descriptor1(i,:),n2,1);
            despxor = xor(descriptor1rep,descriptor2);
            despxor = sum(despxor,2);
%             [~,id] = min(despxor);
%             matching_pairs(i,:) = [i, id];
            [val, id] = sort(despxor,'ascend');
            if val(1) < 0.7 * val(2)
                matching_pairs(i,:) = [i, id(1)];
            else
                matching_pairs(i,:) = [i, inf];
            end
        end
       
    varargout{1} = matching_pairs;
end
