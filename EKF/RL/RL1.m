clear
close all
rng(100)

%  Set to true, to resume training from a saved agent
% resumeTraining = true;
% % Set ResetExperienceBufferBeforeTraining to false to keep experience from the previous session
% agentOpts.ResetExperienceBufferBeforeTraining = ~(resumeTraining);
% if resumeTraining
%     % Load the agent from the previous session
%     sprintf('- Resume training of: %s', PRE_TRAINED_MODEL_FILE);   
%     load(PRE_TRAINED_MODEL_FILE,'saved_agent');
%     agent = saved_agent;
% end
    % Create a fresh new agent
%     agent = rlDDPGAgent(actor, critic, agentOpts);
% end
% Train the agent
% trainingStats = train(agent, env, trainOpts);


record={0,0,0};
save('Qs.mat','record')
%% Environment
% actInfo = rlNumericSpec([1 21]);
actInfo = rlNumericSpec([1 6]);
actInfo.Name = 'actionspace';
actInfo.Description = 'first 6: Qdiag,\ 7-22: Qoffdiag,\ 23-27: P0';
% actInfo.LowerLimit = [zeros(1,6) , -100*ones(1,15) , zeros(1,6)];
% actInfo.LowerLimit = [zeros(1,6) -0.5*ones(1,15)];
% actInfo.UpperLimit = [10 10 1 10 10 1 0.5*ones(1,15)];

% actInfo.LowerLimit = [0.00001*ones(1,6) -0.01*ones(1,15)];
actInfo.LowerLimit = 0.0001;
% actInfo.LowerLimit = -0.8*ones(1,1);
% actInfo.UpperLimit = [2 2 0.2 4 4 0.2 0.01*ones(1,15)];
actInfo.UpperLimit = [0.1 0.1 0.1 0.1 0.1 0.1];
% actInfo.UpperLimit = 0.8*ones(1,1);

obsInfo = rlNumericSpec([1 7000]);
obsInfo.Name = 'observationspace';
obsInfo.Description = '6x1000 EXK; 1x1000 NEES';


% obsPath = [featureInputLayer(7000, 'Normalization','none','Name','myobs') 
%     fullyConnectedLayer(1,'Name','obsout')];
% actPath = [featureInputLayer(21, 'Normalization','none','Name','myact') 
%     fullyConnectedLayer(1,'Name','actout')];
% 
% comPath = [additionLayer(2,'Name', 'add')  fullyConnectedLayer(1, 'Name', 'output')];
% 
% net = addLayers(layerGraph(obsPath),actPath); 
% net = addLayers(net,comPath);
% 
% net = connectLayers(net,'obsout','add/in1');
% net = connectLayers(net,'actout','add/in2');
% 
% % plot(net)
% 
% critic = rlQValueRepresentation(net,obsInfo,actInfo, ...
%     'Observation',{'myobs'},'Action',{'myact'});
% 
% % v = getValue(critic,{rand(7000,1)},{rand(21,1)})
% 
env = rlFunctionEnv(obsInfo,actInfo,'EKFtrain','resetfxn');
% validateEnvironment(env)
%% Q-Learning Agent
% https://www.mathworks.com/help/reinforcement-learning/ug/train-q-learning-agent-to-solve-basic-grid-world.html
% qTable = rlTable(getObservationInfo(env),getActionInfo(env));
% qRepresentation = rlQValueRepresentation(qTable,getObservationInfo(env),getActionInfo(env));
% qRepresentation.Options.LearnRate = 1;

% agentOpts = rlQAgentOptions;
% agentOpts.EpsilonGreedyExploration.Epsilon = .04;
% agent = rlQAgent(critic,agentOpts)
% initOpts = rlAgentInitializationOptions('NumHiddenUnit',128);
% agent = rlSACAgent(obsInfo,actInfo,initOpts);
% qAgent = rlQAgent(critic,agent);


% opt = rlRepresentationOptions('UseDevice',"gpu");
opt = rlDDPGAgentOptions;
opt.NoiseOptions.Variance = 0.09;
opt.NoiseOptions.MeanAttractionConstant = 0.003;
opt.NoiseOptions.VarianceMin = 0.001;
opt.NoiseOptions.VarianceDecayRate = 0;
opt.NoiseOptions.Mean = 0;
% opt.NoiseOptions.VarianceDecayRate = 1e-5;
agent = rlDDPGAgent(obsInfo,actInfo,opt);
% agent = rlSACAgent(obsInfo,actInfo);
%% Train Q-Learning Agent
trainOpts = rlTrainingOptions;
trainOpts.MaxStepsPerEpisode = 1;
trainOpts.MaxEpisodes= 1000;
trainOpts.StopTrainingCriteria = "AverageReward";
trainOpts.StopTrainingValue = 1000000;
trainOpts.ScoreAveragingWindowLength = 30;
trainOpts.SaveAgentCriteria = "EpisodeReward";
trainOpts.SaveAgentValue = 1000;

% trainOpts.UseParallel = true;
% trainOpts.ParallelizationOptions.Mode = "sync";
% trainOpts.ParallelizationOptions.DataToSendFromWorkers = "Experiences";
% trainOpts.ParallelizationOptions.StepsUntilDataIsSent = -1;

% trainingStats.Information.HardwareResource.actorDevice(1) = "gpu";
% trainingStats.Information.HardwareResource.criticDevice(1) = "gpu";
% doTraining = false;
% figure; hold on
trainingStats = train(agent,env,trainOpts);
% hold off






