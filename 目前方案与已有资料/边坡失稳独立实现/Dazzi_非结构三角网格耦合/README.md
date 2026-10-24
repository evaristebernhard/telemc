# Dazzi 非结构三角网格独立实现

本目录用于整理 Dazzi 等（2019）方法，并准备后续 TELEMAC-2D 非结构三角网格边坡失稳实现。

- 方法理解文档：`docs/Dazzi_方法理解与非结构三角网格改造约束.docx`
- 原始论文：仓库根目录的 `Dazzi 等 - 中英对照版.pdf`
- `src/`、`scripts/`、`config/`、`tests/` 用于后续独立实现。
- 现有局部案例、原始几何和基线 Fortran 文件不在本目录内修改。

核心参数必须分开管理：`phi_crit_dry`、`phi_crit_wet`、`phi_deposition`。
