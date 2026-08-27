# Paper 2 D45：endpoint—landing—胃温摄食的前后时间方向检验

日期：2026-08-22  
状态：`OUTCOME_INFORMED_DIRECTIONAL_GAP_COMPLETION__LOCKED_BEFORE_NEW_SIGNED_LAG_RESULTS`

## 1. 触发与纠错

Uesaka旧主结果把`landing in (0,600] s after endpoint`预先定义为结局，只能回答future landing，不能据此回答landing在endpoint前后哪一侧富集。RFBO的`dive_near`则是endpoint `±4 min`时间邻近标签。两者不能直接比较。用户要求对没有计算过的前后方向老实计算，并把已有胃温真值纳入，而不是复述旧结论。

本实验是已知边际结果后的方向缺口补齐，不伪称独立盲发现；所有输入、尺度、时间窗和门在首次生成新的signed-lag统计前冻结。

## 2. 冻结输入

### Uesaka landing

- 10--100 m冻结middle endpoints：`PAPER2_ALBATROSS_ECOLOGY/UESAKA_RD10_100_CHL_LANDING_SENSITIVITY_V1/formal/rd10_100_middle_events_with_gps_outcome.csv.gz`；
- 50--2000 m冻结全尺度事件：`PAPER2_ALBATROSS_ECOLOGY/UESAKA_STEP_METHOD_RD_DELTA_LANDING_SENSITIVITY_V1/formal/rd_cross_delta_events.csv.gz`；正式时使用同目录与其键一一一致的`method_events_with_gps_outcome.csv.gz`读取冻结GPS双侧覆盖字段，只取`phase=middle`；
- 冻结GPS状态转换：同目录`gps_state_transitions_reproduced.csv.gz`中的`landing_gps_conservative`。
- 双侧完整middle覆盖：同目录`phase_runs.csv.gz`；要求endpoint index距其所在middle run两端各至少660个原生1-s位置，并同时通过冻结`primary_risk_eligible_gps`，避免记录/阶段边缘制造前后差。首次validation-only在未计算真实方向效应时发现index推算run时刻跨尺度最大spread为49 s，故在原600点上结果前固定增加60点安全垫，保证实际双侧仍超过600 s；统计窗本身不变。
- 尺度完整保留`10/20/30/40/50/75/100/150/200/300/500/750/1000/1500/2000 m`；10--100 m优先使用第一张结果表，150--2000 m使用第二张，不按方向结果选delta。

### D09胃温与landing

- `PAPER2_ALBATROSS_ECOLOGY/D09_DMS_SOURCE_FOOD_OUTCOME_V1/analysis_formal_v1/primary_landing_sample.parquet`中的冻结2,546个主QC landing及胃温记录覆盖；
- 同目录`ingestion_with_dms.parquet`中的403个作者胃温摄食事件；本实验只读行为时间/ID/覆盖，不读DMS或CHL效应。

### D02胃温与冻结RD

- `codex_runs/d02_eight_str_full_analysis_20260629_0750/.../D02_behavior_episodes_primary.parquet`中的`ingestion_candidate_episode`；
- `codex_runs/albatross_all_datasets_external_validation_audit_20260628_185339/.../D02_boundary_times_long.parquet`中的冻结RD100/RD200 endpoint。
- D02是黑眉信天翁小样本候选胃温事件；它只作真正endpoint—胃温方向桥，不与D09漂泊信天翁作者摄食事件混池。

## 3. Uesaka主检验：endpoint中心的landing前后

对每个middle RD endpoint，在同deployment内计算所有保守landing的有符号lag：

`lag = landing_time - endpoint_time`。

- `lag>0`：landing发生在endpoint之后；
- `lag<0`：landing发生在endpoint之前。

主窗为对称`10 min`，并固定四个互斥lag bins：`0--30 s`、`30--120 s`、`120--300 s`、`300--600 s`，前后镜像。每个endpoint在每个方向/bin只记是否至少一个landing，避免同一侧多次transition重复加权；同时保存全部pair和最近lag。

主统计量：

1. 每尺度`P(post10)-P(pre10)`；
2. post10与pre10各自相对结构null的excess；
3. 各互斥bin的post-minus-pre，定位方向从何时开始。

不使用旧`landing_10m`列作新结果，只用原始冻结事件时间重新构造两侧。

结构null：在每个deployment内把完整landing时间序列共同循环平移一个随机offset，保留landing之间的间隔、endpoint序列、deployment时长和事件数，5,000次；不逐landing独立打乱。对每尺度报告实际、null均值、excess、方向性随机化p。整只bird 20,000次bootstrap给raw post-minus-pre和null-centered direction的95% CI。

方向门：某尺度称`post-dominant`必须`post-minus-pre`的bird CI下界`>0`，null-centered CI下界`>0`，且15尺度Holm后方向p`<=.05`。至少相邻两个尺度通过才称方向尺度族。若post/pre均高于null而差异失败，写“双侧邻近”；若仅post过门，写“endpoint先于landing”；若仅pre过门，写“landing先于endpoint”。

单侧邻近门对15尺度×post/pre共30项作Holm；某侧须null-centered bird CI下界`>0`且Holm `p<=.05`才称该侧高于结构null。方向门与单侧门分别校正，不混为同一家族。

解释限制：endpoint在飞行点上定义，前后机会并非随机试验；即使post-dominant也不能自动称endpoint导致landing。phase-shift只检验时间锁定，不能消除未观测共同状态。

## 4. D09：landing中心的作者胃温摄食方向

对具有胃温记录`±60 min`完整覆盖的主QC landing，计算同deployment作者胃温事件：

`lag = ingestion_time - landing_time`。

固定镜像bins为`0--5 min`、`5--15 min`、`15--60 min`，主统计`P(post60)-P(pre60)`。用deployment内胃温事件序列共同循环平移5,000次，整只bird bootstrap 20,000次。单一主门要求raw及null-centered post-minus-pre CI下界`>0`且随机化`p<=.05`。

另直接连接ingestion是否落在该landing开启的wet bout `[bout_start,bout_end]`，报告landing→ingestion lag分布；这回答真实摄食相对海面转换的顺序，不把没有记录摄食的landing称失败。

## 5. D02：冻结RD endpoint中心的胃温候选方向

仅用RD100、RD200 endpoint与`ingestion_candidate_episode` onset；同bird/trip内计算：

`lag = ingestion_time - endpoint_time`。

镜像主窗15 min，互斥bins`0--5`、`5--15 min`。事件须位于该trip endpoint时间范围内且两侧各有15 min覆盖。对每种RD报告endpoint后摄食与endpoint前摄食的事件比例、bird bootstrap CI及trip内共同循环相位null。RD100/RD200两项Holm；因胃温事件少，覆盖不足按`not_estimable`处理，不用D09救回。

## 6. 共同边界

- D09的胃温摄食与Uesaka不是同一deployment，不能给Uesaka landing逐次贴成功标签；
- D02为另一物种、候选事件和小样本，只是方向复现层；
- 阳性不证明CHL、视觉、嗅觉或猎物cue触发endpoint；
- 不重跑RD、CHL、阶段、landing识别、胃温识别或Lévy分类；不修改Paper 1。

## 7. 执行约束

- validation先只核对键、时区、唯一性、尺度、记录覆盖和合成时滞，不输出真实方向效应；
- 正式CPU7、BLAS单线程、CPU1禁止；20,000 bird bootstrap、5,000 phase shifts；
- 无下载、无网络、无外部联系；结果先在对话中报告，再同步Paper 2文档。

冻结输入SHA-256：

- Uesaka 10--100：`86fcaa77f31e2ebc1489cc78e011a42eca7b86704bd2a7725140a1c26cbb02c9`
- Uesaka 50--2000：`848c479931cb198b03998d02681d4382c576e3853bc91187a1665b2946108ac1`
- Uesaka 50--2000 GPS outcome：`6b6b66fe931c11f093233fb522afeb47140b3d8c13ce3038cdb090285a852ee0`
- Uesaka landings：`1e3efb631cc4dd3393217f246e3ad85ebfa183b1cf819700950e8ef43acf7cdb`
- Uesaka phase runs：`46fce2ab1791b533f33b4d6e1d2787f023f60da25471a15f05ca07a287e1b90c`
- D09 landings：`1167d104e43e28dd8905f1cd73754d67778e0f17720d295d1608a24747897340`
- D09 ingestions：`9560190b566fcbb80513393f88b50bc07b1c4c1f92e6ed4f98b3f70e708b8a14`
- D02 episodes：`ddc0ddf160bbbc12ae9220a67241f0c14890406a4024b2d71a04ab94d20e59f9`
- D02 boundaries：`0c114b2849689f7381a85ed62a87034fa332436e896cfde4b0a8105355de8e73`

## 8. 2026-08-22 正式v1后的技术纠错（新D09结果前冻结）

`formal_v1`完成后，D09出现“全部±60 min均为0”，与原表中肉眼可见的同deployment相邻秒级事件矛盾。只读核查确认：

- landing时间列为`datetime64[ns, UTC]`；
- ingestion时间列为`datetime64[us, UTC]`；
- 旧实现直接取各自整数再统一除以`1e9`，把ingestion epoch缩小了1,000倍；
- 这是代理实现错误，不是生物学零结果。

处理决定：

1. 所有datetime先显式归一到`datetime64[ns, UTC]`再转epoch秒；
2. `formal_v1`的Uesaka结果有效（输入时间均为ns），D02输入时间也均为ns；保留它们，不重跑；
3. 只重算受影响的D09 landing↔作者胃温方向，仍使用本文件第4节已冻结的窗、bins、5,000 phase shifts、20,000 bird bootstrap、固定种子与判据；
4. 旧D09零结果标记为`TECHNICALLY_INVALID_DATETIME_UNIT_MISMATCH`，禁止引用。
