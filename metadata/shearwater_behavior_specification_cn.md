# Paper 2 D107：短尾鹱 last-passage CHL 足迹与未来捕食—休息状态直接连接

日期：2026-08-25  
状态：`LOCKED_BEFORE_DOWNLOADING_OR_READING_ANY_SHORT_TAILED_SHEARWATER_CHL_VALUE`

## 1. 唯一问题

D104证明了若干鸟类系统的step-relative CHL端点富集主要位于strict radial records中的最后记录选择：

`L = I(last radial record in tail) - mean(all radial-record tail flags)`。

D106在RFBO中发现潜水尝试事件的`L_union`在1--2 km为正，但它不强于只有落水而没有潜水的事件，因此不能称捕食特异。D107使用一套独立、公开、带逐秒加速度行为标签的短尾鹱数据，检验：**同一冻结RD事件的未来5 min越偏向捕食而非休息，其last-passage CHL双尾选择是否越强。**

阳性只允许解释为“last-passage环境足迹与捕食状态同事件相连”；不能解释为鸟直接感知CHL、CHL触发转弯、捕获成功、视觉或嗅觉机制已经识别。

## 2. 结果盲规划及冻结输入

冻结本文件前没有下载、打开或计算短尾鹱2012年CHL值。只完成了运动几何、公开行为覆盖和最小下载清单规划：

- 公开行为原始文件SHA256：`09ee71be7a0be335ecbf5b0626ca4339bc7f0daa4c50af2c8a62af79860b0964`；
- 已完成实验冻结的phase-0、5-min支持、RD-P端点表SHA256：`40f0bea29d9f868121641a3270a78d8864390b3e35da7753698091278c974b43`；
- 结果盲事件几何目录：`PAPER2_ALBATROSS_ECOLOGY/D107_SHORT_TAILED_LAST_PASSAGE_CHL_FORAGING_REST_V1/plan_v1/geometry_event_catalog.csv.gz`，SHA256 `d897562ef78d5dd77161dfe85a8ce3aa8245386dcdb6605a8c4074ae2a86ba9a`；
- 下载清单SHA256：`e051affbd7c6d0136a04e773e18c7756bdb4b4cf0fc8038725f06c08e3d69a20`；
- 产品：Copernicus Marine `cmems_obs-oc_glo_bgc-plankton_my_l4-gapfree-multi-4km_P1D`，日尺度、原生约1/24度、变量`CHL`；
- 只下载10个观测日期、24个日期×1度tile；不得因结果方向补日期、补区域或更换产品。

不得重跑RD，不得使用逐秒插值坐标生成新端点。运动只使用既有phase-0规则5-min支持；公开逐秒坐标只用于在相邻硬件GPS支持间汇总既有行为窗，不冒称原生1-Hz GPS。

## 3. 冻结尺度、事件与环境序列

固定`delta={500,1000,2000,5000} m`。10 km已在行为桥实验中覆盖不足，不得结果后加入；10--100 m不受5-min GPS支持。

每个run内按冻结端点递归切分：第一步起点为该run索引0，之后一步起点为上一个冻结端点，终点为当前冻结端点。对每个事件：

1. 用每个5-min支持点的日期与坐标匹配最近的原生CHL网格中心；不做双线性插值或伪造亚网格分辨率；
2. 同日期、同原生网格的连续支持点折叠为一个cell-run；以`log(CHL)`为环境值；
3. 事件必须具有至少3个有限且非恒定cell-runs；
4. 从事件起点计算径向距离，冻结strict radial records；记录点映射到至少2个不同cell-runs；最后一个record必须为冻结RD endpoint；
5. 缺值、日期不覆盖、几何恒等式失败均排除并逐事件记录原因，不允许插补CHL。

冻结规划中共有314个端点、254个仅按几何可用的事件；四尺度几何可用事件为78/73/65/38，均覆盖9只个体。CHL读取后的最终覆盖允许因缺值下降，但不得因效应方向改变排除规则。

## 4. 冻结行为量与last-passage量

行为完全复用既有短尾鹱行为桥定义。每个冻结endpoint未来`[0,300) s`：

- `p_forage`：`diving OR sforaging`秒数/300；
- `p_rest`：`resting`秒数/300；
- `specificity = p_forage - p_rest`。

`diving/sforaging`表示加速度识别的捕食尝试/觅食行为，不表示捕获成功。

CHL在每个事件cell-run序列内以平均秩定义`high: rank>0.8`、`low: rank<=0.2`、`union=high OR low`。主环境量是：

`L_union = I(endpoint cell-run in union) - mean(radial-record cell-run union flags)`。

并完整保存`E_high/E_low/E_union`、`R_high/R_low/R_union`、`L_high/L_low/L_union`，验证逐事件`E=R+L`。

## 5. 主估计量

每尺度先在每只生物个体内中心化行为特异性和`L_union`。主交互量为个体内线性斜率：

`beta_i = cov(L_union, specificity) / var(specificity)`，

再对符合条件的个体等权平均。个体需至少3个最终可用事件且`specificity`有非零方差；尺度覆盖门为至少6只个体、至少30个事件。

为避免“正斜率只是休息事件更负、但捕食事件本身没有last-passage选择”，第二主量是`forage-dominant`事件（`specificity>0`）中先对个体内事件等权、再对个体等权的`L_union`均值。该量覆盖门为至少5只个体、15个事件。

分类的鸟内`forage-dominant minus rest-dominant`差为次级估计；由于5000 m在结果盲规划中只有4只配对鸟，不承担主门。

## 6. 共同相位null与置信区间

零模型必须保留真实RD端点、事件边界、strict-record位置、行为值、个体结构、CHL值分布、自相关和趋势：

- 每个`scale × support run`把其递归事件cell-run序列按事件顺序连接成master token序列；
- 每次置换对同一support run内全部事件施加同一个非零circular CHL phase，再按原事件宽度和record位置重算全部`L`；
- 行为标签不移动；同一事件的运动几何不移动；
- 对每尺度999次有效共同phase，固定正式seed `1070825`；不得删除不利phase；
- 主斜率和forage-dominant组内量各在四尺度内独立做单侧phase p并Holm校正。

观察量的95% CI以生物个体为cluster、20,000次等权bootstrap；不把deployment或事件当独立个体。正式结果必须报告实际事件数、个体数、run数、有效phase数、重复tile值一致性、全部输入/脚本/输出哈希和CPU affinity。

## 7. 冻结裁决门

单尺度主斜率通过要求：

1. 覆盖门通过；
2. `mean(beta_i)>0`；
3. individual-bootstrap 95% CI下界`>0`；
4. 四尺度Holm校正phase `p<=.05`。

同尺度forage-dominant `L_union`通过要求：覆盖通过、观察值`>0`、individual-bootstrap 95% CI下界`>0`、四尺度Holm phase `p<=.05`。

只有至少两个相邻尺度同时通过上述两个门，才裁决`DIRECT_FORAGING_STATE_LINK_SUPPORTED`。否则裁决`NO_DIRECT_FORAGING_STATE_LINK`，不得改变行为窗、tail cut、网格提取、尺度、个体权重、phase方式或裁决门救回。

## 8. 次级分解和解释边界

必须完整报告但不能救主结果：

- `L_high`与`L_low`的连续斜率及forage-dominant组内量，按`2 tails × 4 scales`校正；
- `E`与`R`的union斜率，判断关联是否特异于last-passage而非普通endpoint excess或record sampling；
- `p_forage`和`p_rest`分开替代specificity；
- 过去5 min同定义行为量以及future-minus-past方向，仅用于区分持续觅食状态与单向行为起点；
- 分类的forage-dominant、rest-dominant和两者鸟内差。

主门通过仍不得写：CHL是鸟直接感知的cue、CHL因果触发转弯、识别了视觉/嗅觉、证明捕获成功、或已单独生成/推翻Lévy。主门失败则说明当前两套同事件行为数据都没有把总体last-passage CHL足迹确定为捕食特异机制，必须继续公开同步行为数据，而不能用跨文献拼接宣告机制完成。

正式计算固定非CPU1单核心、BLAS单线程；只用公开数据，不外部联系，不修改Paper 1。
