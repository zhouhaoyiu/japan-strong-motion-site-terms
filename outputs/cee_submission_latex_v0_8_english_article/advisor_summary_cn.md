# 给导师的工作说明

## 题目

Long-period site residuals and response-spectrum hazard sensitivity in Japanese strong-motion records

## 核心问题

本文研究日本强震记录中 1--3 s 长周期场地残差是否具有稳定的非遍历特征，以及这种台站项是否会改变官方 J-SHIS 长周期响应谱坐标。

## 已完成工作

1. 建立了 J-SHIS/NIED 强震记录残差数据集，使用 333,808 条 K-NET 和 KiK-net 记录计算 MF2013 地震动模型残差。
2. 证明 MF2013 的 D1400/AVS30 场地项主要改进 1--3 s 周期段：PGA、SA(0.3 s)、SA(1.0 s)、SA(3.0 s) 的平均绝对残差降低分别为 1.6%、3.9%、16.1% 和 21.2%。
3. 证明 SA(3.0 s) 台站残差与 D1400 的相关性在官方场地校正后显著降低，Spearman 相关系数由 0.614 降至 -0.025，说明官方场地项吸收了主要长周期盆地深度结构。
4. 建立台站留出和空间分块验证模型，SA(1.0 s) 和 SA(3.0 s) 的空间分块 RMSE 分别降低 38.9% 和 44.3%。
5. 增加可解释场地变量模型，不使用台站编号，仅使用 D1400、Dbase、AVS30、VS20、钻孔、地形、火山前缘距离和区域项。该模型对 SA(3.0 s) 台站项的 RMSE 相对零校正降低 34.7%。
6. 将 SA(3.0 s) 台站校正因子施加到官方 J-SHIS 响应谱网格。50 年 10% 超越概率下，匹配台站网格的 SA(3.0 s) 中位数由 0.076 g 变为 0.048 g。
7. 完成 10,467 个网格的全周期响应谱传播审计，50 年 10% 水平下 0.3 s、1.0 s、3.0--5.0 s 的中位倍率分别为 1.045、0.834 和 0.658。
8. 完成 J-SHIS 公开源参数和采样网格一致性检查：421 个官方 PRM CSV、200 个 shapefile 图层完成清点，10,467 个抽样网格和 41,868 条官方阈值行通过检查。

## 文章主线

日本强震记录中的 1--3 s 非遍历场地残差可以被公开场地变量预测，并且足以改变官方长周期响应谱坐标。文章不声称完整复现 J-SHIS 官方生产级 PSHA；官方空间相关、逻辑树聚合、模型不确定性和生产系统实现仍保留在 J-SHIS 官方产品中。

## 投稿定位

建议目标期刊为 Communications Earth & Environment。文章的强项是公开数据、强震记录规模、台站留出验证、官方响应谱坐标影响和可追踪的公开参数检查。主要风险是审稿人是否接受“官方响应谱敏感性 + 公开参数一致性检查”作为足够强的危险性证据。

## 当前材料

- 英文主文：`main.pdf`
- 补充材料：`supplementary_information.pdf`
- Cover letter：`cover_letter_cee.md`
- 投稿检查清单：`submission_closure_checklist_cee.md`
- 代码和派生表仓库：`https://github.com/zhouhaoyiu/japan-strong-motion-site-terms`
