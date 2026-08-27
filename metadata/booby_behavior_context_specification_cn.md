# Paper 2 D106：RFBO last-passage 分量与同事件真实捕食意图直接连接

日期：2026-08-25  
状态：`LOCKED_BEFORE_ANY_D104_L_BY_BEHAVIOR_EFFECT_WAS_CALCULATED`

## 1. 唯一问题

D104已经把跨系统的step-relative CHL双尾定位到未来drawdown从strict radial records中选择出的最后记录点：

`L = I(last radial record in tail) - mean(all radial-record tail flags)`。

它回答了endpoint为什么在数学上特殊，却没有回答这种last-passage选择是否与真实觅食行为有关。D106只检验：**RFBO同一RD事件中的`L`是否在真实TDR潜水尝试附近增强，并且强于只有落水而没有潜水尝试的事件。**

若不强于仅落水，不能把总体`L`解释成捕食特异机制；若增强，也只能说明last-passage环境足迹与捕食意图同事件相连，不能证明CHL是鸟的感觉cue或因果刺激。

## 2. 已知与未知，防止结果后叙事

冻结本文件前已经知道：

- D104中RFBO总体`L_union`在250/500/1000/2000 m均通过强phase null；
- D28的普通endpoint excess `E`在`dive_near`内部多尺度为正，但`dive_near-wet_only_near`的union/low没有形成稳定直接差，500 m high反而较低；
- 潜水发生在endpoint前后两侧，方向不对称失败，因此endpoint只能视为持续局地觅食状态中的renewal标记，不能先验写成单向刺激起点；
- D106覆盖审计只读取键和标签，未读取任何`L/R`值：四尺度分别有88/82/72/61只鸟同时贡献`dive_near`和`wet_only_near`事件。

冻结前**没有**计算、查看或推断任何`L_high/L_low/L_union × behavior_class`组均值、行为差、phase p或置信区间。

## 3. 冻结输入与标签

- D104 Family A逐事件表：`PAPER2_ALBATROSS_ECOLOGY/D104_LAST_PASSAGE_CHL_DECOMPOSITION_V1/formal_v2_phase999/family_a_event_metrics.csv.gz`，SHA256 `aeb4cfb5b3fb8a5e38524e4a7f597a0e966fc9cc38a440591d987d8d0c86e23b`；
- D28冻结行为标签：`PAPER2_ALBATROSS_ECOLOGY/D28_RFBO_RD_ENDPOINT_STEP_CHL_TWO_TAIL_V1/formal_v1/event_two_tail_metrics.csv.gz`，SHA256 `017fceace0c54e57d315dd1d7fb96bd6f21ee3d6423a73bd6ba9248a8c92d9c4`；
- D104原始RFBO GPS＋CHL、movement-only RD事件和deployment→BandNum映射只用于重建phase null，不重跑、重选或修改RD；
- 固定尺度为250/500/1000/2000 m，4000 m不在D104 Family A冻结尺度族，不得结果后加入；
- `dive_near`：endpoint前后`±240 s`内有GPS-linked、校正深度`>=1 m`的TDR潜水尝试；`wet_only_near`：同窗口有wet-no-dive且无潜水；`no_tdr_near`：同窗口无同步TDR行为事件。标签优先级与D28完全一致。

潜水代表觅食/捕食尝试，不要求捕获成功；wet-only代表坐水/湿水而无已检测潜水，不等于生理上确认休息。

## 4. 主估计量

每个尺度只用D104 Family A可分解且与D28一一连接的RFBO事件。先在每只`BandNum`内部对事件等权，再对鸟等权。

两个并列且缺一不可的主量都只看`L_union`：

1. `DIVE_L_UNION`：`dive_near`事件中鸟等权的`L_union`均值；
2. `DIVE_MINUS_WET_L_UNION`：同时具有两类事件的鸟内，`mean(L_union|dive_near)-mean(L_union|wet_only_near)`，再对配对鸟等权。

第1项回答捕食尝试附近是否存在正last-record selection；第2项回答它是否超出普通落水背景。不能用第1项阳性替代第2项。

覆盖门：每尺度`dive_near`与`wet_only_near`各至少100 events、30 birds，且至少30只配对鸟。所有覆盖判定只按键与标签，不看效应方向。

## 5. 结构null与不确定性

- 完全复用D104的Family A cell-run序列、strict radial-record位置和每segment共同的非零circular CHL phase；同一segment内所有事件共享相位，保留运动几何、事件边界、record位置、CHL自相关、趋势、窗口长度、行为标签位置及同segment依赖。
- 每个尺度999次有效phase，固定seed `1060825`；不得删掉不利phase。
- 观察效应以BandNum为单位bootstrap 5,000次；组内量重采有该类事件的鸟，行为差重采同时有两类事件的配对鸟，报告双侧95% CI。
- 对`DIVE_L_UNION`和`DIVE_MINUS_WET_L_UNION`分别计算观察量相对phase-null的单侧p；每个量的四尺度各自构成一个Holm family。两组family不得合并后挑较小者，也不得用未校正p裁决。

## 6. 冻结裁决门

单尺度`DIVE_L_UNION`通过要求：覆盖通过、观察值`>0`、bird-bootstrap CI下界`>0`、四尺度Holm phase `p<=.05`。

单尺度`DIVE_MINUS_WET_L_UNION`通过要求：覆盖通过、观察差`>0`、配对bird-bootstrap CI下界`>0`、四尺度Holm phase `p<=.05`。

只有同一组至少两个相邻尺度同时通过上述两个门，才裁决`DIRECT_FEEDING_INTENT_LINK_SUPPORTED`。否则裁决`NO_DIRECT_FEEDING_INTENT_SPECIFICITY`，不得更换窗口、深度阈值、尺度、tail cut、个体权重或行为对照救回。

## 7. 必须完整报告但不能救主结果的分解

- `L_high`、`L_low`分别重复组内与`dive-wet`估计，各自按`2 tails × 4 scales`校正；
- `dive_near-no_tdr_near`为次级对照，三尾按`3 × 4`校正；
- `R_high/R_low/R_union`及`E_high/E_low/E_union`重复行为差，用于判断若有差异是否特异于last-record selection，而非一般radial-record sampling或原始endpoint excess；
- 报告每尺度/行为的events、birds、segments和配对鸟数，以及精确一一连接、`E=R+L`恒等式、CPU affinity、输入与脚本哈希。

任何单尾、单尺度或未校正方向都只能标为描述性，不得升级为机制结论。

## 8. 允许与禁止解释

若主门通过，允许写：RFBO中与真实潜水尝试同事件的last-passage CHL足迹强于仅落水事件，支持环境条件化的觅食状态边界。仍不得写CHL被直接感知、CHL触发转弯、潜水成功捕获、视觉/嗅觉已区分、或已生成/推翻Lévy机制。

若主门失败，D104仍是稳健的路径—环境联合结构，但当前RFBO同事件行为不支持把它解释成捕食特异机制；必须继续寻找另一套公开、同步行为真值，而不是以跨数据文献拼接宣布完成。

正式计算固定单一非CPU1核心、BLAS单线程；不修改Paper 1，不外部联系，只用已公开数据。
