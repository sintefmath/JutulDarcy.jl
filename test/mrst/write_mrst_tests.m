%% MRST (www.sintef.no/mrst) script to generate Jutul input test files
% More user friendly scripts are available with MRST in the jutul module.
mrstModule reset ad-core ad-blackoil deckformat jutul test-suite
cases = {'spe1', 'spe3', 'spe9', 'egg'};
n = numel(cases);
paths = cell(n, 1);
for i = 1:numel(cases)
    name = cases{i};
    switch name
        case 'spe1'
            pth = getDatasetPath('spe1');
            deck_path  = fullfile(pth, 'BENCH_SPE1.DATA');
        case 'spe3'
            pth = getDatasetPath('spe3');
            deck_path  = fullfile(pth, 'BENCH_SPE3.DATA');
        case 'spe9'
            pth = getDatasetPath('spe9');
            deck_path  = fullfile(pth, 'BENCH_SPE9.DATA');
        case 'egg'
            deck_path = getDeckEGG();
        otherwise
            error('No such case.')
    end
    [state0, model, schedule, nls] = initEclipseProblemAD(deck_path, 'ReorderStrategy', reorder);
    [ws, states] = simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls);
    ref.bhp = getWellOutput(ws, 'bhp');
    if model.water
        ref.qWs = getWellOutput(ws, 'qWs');
    end
    if model.oil
        ref.qOs = getWellOutput(ws, 'qOs');
    end
    if model.gas
        ref.qGs = getWellOutput(ws, 'qGs');
    end
    ref.names = {schedule.control(1).W.name};
    paths{i} = writeJutulInput(state0, model, schedule, name, 'extra', struct('mrst_solution', ref));
end
