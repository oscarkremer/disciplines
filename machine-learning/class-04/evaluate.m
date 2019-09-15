function classification = evaluate(result, testY) 
    classification = zeros(2,2);  
    for i=1:500
        if result(i) == -1 && testY(i) == -1
            classification(1,1) = classification(1,1) + 1;
        end
        if result(i) == -1 && testY(i) == 1
            classification(1,2) = classification(1,2) + 1;
        end
        if result(i) == 1 && testY(i) == -1
            classification(2,1) = classification(2,1) + 1;
        end
        if result(i) == 1 && testY(i) == 1
            classification(2,2) = classification(2,2) + 1;
        end
    end
end
