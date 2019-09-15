function probability = prob(word, num) 
    length = size(word);length=length(1);
    probability = 1;
    for i=1:length
        probability = probability*((word(i))^num(i));
    end
end
