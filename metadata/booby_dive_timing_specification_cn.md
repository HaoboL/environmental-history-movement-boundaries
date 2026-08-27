# D126 RFBO last-passage × 潜水前后方向：预注册

日期：2026-08-26

## 新问题与既有披露

D28-S2C/D49已经证明TDR潜水在RD endpoint前后两侧均富集，前后不对称性失败；D49的原始endpoint CHL high/low也不按pre/post分化。D104后来才把endpoint excess分为record-frontier `R`和last-passage `L`，既有时间方向分析没有检验`L`。

D126只回答一个未重复问题：`L_union`在endpoint之后将潜水（`post_only`）时是否强于潜水后才出现endpoint（`pre_only`）。这区分“捕食尝试前的接近/状态边界”与“捕食尝试后的路径重置”，不重跑RD、TDR、CHL或既有D49。

因D49原始CHL方向结果已披露，D126属于D104之后的机制分解，不冒充完全盲法发现。

## 冻结输入

- D104 formal_v2 Family A RFBO冻结事件，尺度250、500、1000、2000 m。
- D49 formal_v3_final 的RFBO行为方向标签，窗口固定为endpoint前后各240 s，潜水阈值固定为校正深度`>=1 m`。
- `post_only`：窗口内只在endpoint之后有潜水；`pre_only`：只在之前有；`both`和`none_observed`不进入主对比。
- 个体单位固定为D104/D28的`BandNum` biological bird；先鸟内平均，再鸟等权。

## 结果盲覆盖门

资格阶段只读两表的事件键、尺度、鸟和D49 group，不读任何`E/R/L`列。单尺度需：

1. pre和post各至少30个D104合格事件；
2. pre和post各至少20只鸟；
3. 至少15只鸟同时有pre和post事件。

至少两个相邻尺度通过才运行正式效应；失败则停止且不放宽240-s窗口、潜水深度、尺度或门槛。

## 正式主量（覆盖通过后）

`Delta_L_union = mean_bird(L_union | post_only) - mean_bird(L_union | pre_only)`

- 正值：更符合last-passage足迹在潜水前增强；
- 负值：更符合潜水后的路径重置；
- 零：更符合持续foraging/surface state中的双侧renewal，而非固定单向链。

结构null沿用D106：每一segment内对共同CHL token序列作circular phase，所有事件共同移动；999 phase。鸟内配对bootstrap 20,000次，双侧phase p，四尺度Holm。单尺度通过需覆盖、CI排除0且Holm`<=.05`；至少两个相邻尺度同号通过才称时间方向特异。

冻结诊断：`L_high/L_low`、`E_union/R_union`及post/none、pre/none，只用于定位主结果，不能替换主门。

## 解释边界

潜水是尝试而非捕获成功；post/pre是240-s邻近关系而非首次感觉时刻。无论结果如何，都不能证明CHL被感知、嗅觉/视觉、猎物存在、摄食成功或Lévy生成。

