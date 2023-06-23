function stop = saveEvalCount(x, optimValues, state)
    global fvalData;
    fvalData = [fvalData; optimValues.fval];
    stop = false;
end