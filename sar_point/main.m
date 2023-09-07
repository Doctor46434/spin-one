%% 主文件，完成原始信号采集、距离向脉冲压缩、方位向脉冲压缩


[ReWave] = OriginOnePoint;   %% 仿真回波信号

[CompWave] = PulseComp(ReWave);   %%对回波信号进行距离向压缩

[CompWave2] = OrieComp(CompWave);   %%对回波信号进行方向向压缩

% [ReWave] = OriginThreePoint;
% 
% [CompWave] = PulseComp(ReWave);
% 
% [CompWave2] = OrieComp(CompWave);



