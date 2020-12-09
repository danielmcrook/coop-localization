rng(100)

%% Environment
actInfo = rlNumericSpec([1 21]);
actInfo.Name = 'actionspace';
actInfo.Description = 'first 6: Qdiag,\ 7-22: Qoffdiag,\ 23-27: P0';
% actInfo.LowerLimit = [zeros(1,6) , -100*ones(1,15) , zeros(1,6)];
actInfo.LowerLimit = 0;
actInfo.UpperLimit = 100;

obsInfo = rlNumericSpec([7 1000]);
obsInfo.Name = 'observationspace';
obsInfo.Description = '6x1000 EXK; 1x1000 NEES';


obsPath = [featureInputLayer(21, 'Normalization','none','Name','myobs') 
    fullyConnectedLayer(1,'Name','obsout')];
actPath = [featureInputLayer(7000, 'Normalization','none','Name','myact') 
    fullyConnectedLayer(1,'Name','actout')];

comPath = [additionLayer(2,'Name', 'add')  fullyConnectedLayer(1, 'Name', 'output')];

net = addLayers(layerGraph(obsPath),actPath); 
net = addLayers(net,comPath);

net = connectLayers(net,'obsout','add/in1');
net = connectLayers(net,'actout','add/in2');

% plot network
plot(net)



env = rlFunctionEnv(obsInfo,actInfo,'EKFtrain','resetfxn');
validateEnvironment(env)
%% Q-Learning Agent
% https://www.mathworks.com/help/reinforcement-learning/ug/train-q-learning-agent-to-solve-basic-grid-world.html
qTable = rlTable(getObservationInfo(env),getActionInfo(env));
qRepresentation = rlQValueRepresentation(qTable,getObservationInfo(env),getActionInfo(env));
qRepresentation.Options.LearnRate = 1;

agentOpts = rlQAgentOptions;
agentOpts.EpsilonGreedyExploration.Epsilon = .04;
qAgent = rlQAgent(qRepresentation,agentOpts);


%% Train Q-Learning Agent
trainOpts = rlTrainingOptions;
trainOpts.MaxStepsPerEpisode = 50;
trainOpts.MaxEpisodes= 200;
trainOpts.StopTrainingCriteria = "AverageReward";
trainOpts.StopTrainingValue = 11;
trainOpts.ScoreAveragingWindowLength = 30;


% doTraining = false;
trainingStats = train(qAgent,env,trainOpts);









