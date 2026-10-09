# SQP 基础：从约束线性化到 KKT–Newton

> **核心结论**：SQP 每轮用一个二次规划产生候选方向。它同时处理目标下降与约束的一阶修复；线性化可行并不意味着真实非线性可行。

本文是四篇系列的第一篇。先理解方向如何产生，再进入第二篇的二阶结构、第三篇的 Filter 和第四篇的 Funnel。公式统一采用列向量梯度、按行排列的约束 Jacobian，以及 $g(x)\le0$ 的不等式约定。

## 1. SQP 解决什么问题？

考虑光滑非线性规划（NLP）：

$$
\begin{aligned}
\min_{x\in\mathbb R^n}\quad &f(x),\\
\text{s.t.}\quad &c(x)=0,\quad g(x)\le0.
\end{aligned}
$$

其中 $c:\mathbb R^n\to\mathbb R^m$，$g:\mathbb R^n\to\mathbb R^p$。负梯度只告诉我们如何降低目标，不知道哪些方向会破坏约束；可行域通常还是弯曲的。

SQP 的策略是：在当前点把约束线性化，把优化结构近似为二次模型，求出方向后再决定走多远。

$$
\text{NLP}\longrightarrow \text{QP}_0\longrightarrow \text{QP}_1\longrightarrow\cdots,
\qquad x_{k+1}=x_k+\alpha_kd_k.
$$

“序列”意味着在新点重建模型，不是反复求解一个固定 QP。

## 2. Jacobian：变量移动如何改变约束？

把各约束梯度转置后堆成行：

$$
A_k=J_c(x_k)=\begin{bmatrix}\nabla c_1(x_k)^T\\\vdots\\\nabla c_m(x_k)^T\end{bmatrix}\in\mathbb R^{m\times n}.
$$

一阶展开为

$$
c(x_k+d)=c_k+A_kd+O(\|d\|^2),\qquad c_k=c(x_k).
$$

$c_k$ 是当前残差，$A_kd$ 是沿方向移动时约束的一阶变化。于是 $c_k+A_kd=0$ 表示：**局部线性模型预测下一点可行**。

若 $c_k=0$ 且 $A_k$ 满行秩，$A_kd=0$ 描述可行流形的切空间；若 $c_k\ne0$，则 $A_kd=-c_k$ 是修复残差的仿射条件。

> **容易混淆的地方**：可行点的切向条件是 $Ad=0$，不可行点的 SQP 约束通常是 $Ad=-c$。后者不是经过原点的向量空间。

## 3. SQP 的 QP 子问题

记 $q_k=\nabla f(x_k)$、$C_k=J_g(x_k)$，避免把目标梯度与不等式函数 $g$ 混用。

$$
\begin{aligned}
\min_d\quad &q_k^Td+\tfrac12d^TB_kd,\\
\text{s.t.}\quad &c_k+A_kd=0,\\
&g_k+C_kd\le0.
\end{aligned}
$$

常数 $f_k$ 不影响最小解，故省略。理论上的 exact SQP 取

$$
B_k=\nabla^2_{xx}L(x_k,\lambda_k,\mu_k),\quad
L=f+\lambda^Tc+\mu^Tg,\quad \mu\ge0.
$$

它不是只对目标做 Taylor 展开后得到的普通二次模型：$B_k$ 含约束曲率，第二篇会推导其来源。QP 的约束仍然是线性的。

## 4. 单位圆上的第一步

考虑

$$
f(x,y)=(x-2)^2+(y-1)^2,\qquad c(x,y)=x^2+y^2-1=0.
$$

几何上是在单位圆寻找最靠近 $(2,1)$ 的点，真实解为

$$
(x_*,y_*)=\left(\frac2{\sqrt5},\frac1{\sqrt5}\right).
$$

从 $(1,0)$ 出发，取初始乘子 $\lambda_0=0$，于是 $B_0=2I$。当前信息为

$$
q_0=\begin{bmatrix}-2\\-2\end{bmatrix},\qquad
A_0=\begin{bmatrix}2&0\end{bmatrix},\qquad c_0=0.
$$

QP 为

$$
\min_{u,v}\ -2u-2v+u^2+v^2,\qquad 2u=0.
$$

所以 $d_0=(0,1)^T$。它确实沿圆在 $(1,0)$ 的切线移动，但完整步到达 $(1,1)$，真实残差为 $1$。

对任意步长 $\alpha$，有

$$
c(1,\alpha)=\alpha^2,\qquad
f(1,\alpha)=2-2\alpha+\alpha^2.
$$

因此当 $0<\alpha<2$ 时目标下降，约束却出现二阶误差。这不是 QP 算错，而是模型只满足一阶约束。第三篇的接受规则和 SOC 正是为这种差异服务。

## 5. QP 的 KKT 条件

QP 解返回 $d_k$ 和乘子估计 $\widehat\lambda_k,\widehat\mu_k$：

$$
\begin{aligned}
q_k+B_kd_k+A_k^T\widehat\lambda_k+C_k^T\widehat\mu_k&=0,\\
c_k+A_kd_k&=0,\\
g_k+C_kd_k&\le0,\\
\widehat\mu_k&\ge0,\\
\widehat\mu_{k,j}(g_k+C_kd_k)_j&=0.
\end{aligned}
$$

等式情形可以写成 saddle-point 系统：

$$
\begin{bmatrix}B_k&A_k^T\\A_k&0\end{bmatrix}
\begin{bmatrix}d_k\\\widehat\lambda_k\end{bmatrix}
=-\begin{bmatrix}q_k\\c_k\end{bmatrix}.
$$

即使 $B_k\succ0$，整个 KKT 块矩阵通常也是不定矩阵。若 QP 凸，KKT 在适当条件下刻画最优解；若 QP 非凸，仅满足 KKT 不能保证取到最小点。

## 6. SQP 为什么等价于 KKT–Newton？

只考虑等式。原 NLP 的一阶方程为

$$
F(x,\lambda)=\begin{bmatrix}\nabla_xL(x,\lambda)\\c(x)\end{bmatrix}=0.
$$

Newton 线性化得到

$$
\begin{bmatrix}\nabla^2_{xx}L_k&A_k^T\\A_k&0\end{bmatrix}
\begin{bmatrix}d_k\\\Delta\lambda_k\end{bmatrix}
=-\begin{bmatrix}\nabla_xL_k\\c_k\end{bmatrix}.
$$

令 $\widehat\lambda_k=\lambda_k+\Delta\lambda_k$，第一行变成

$$
q_k+B_kd_k+A_k^T\widehat\lambda_k=0.
$$

这正是等式 QP 的 stationarity。因此 exact SQP 的方向可理解为对原 KKT 方程做 Newton。对不等式，局部固定活跃集后可作类似解释，还须处理活跃集改变和互补性。

> **核心结论**：QP 乘子是更新后的乘子估计；Newton 系统里的 $\Delta\lambda$ 是增量。两者相差当前乘子，不能把两个记号混为一谈。

## 7. Exact 与 quasi-Newton / BFGS

Exact SQP 显式计算 Lagrangian Hessian。Quasi-Newton SQP 用梯度差构造近似，以降低二阶导成本。

定义实际接受的位移 $s_k=x_{k+1}-x_k$，使用同一组新乘子构造

$$
y_k=\nabla_xL(x_{k+1},\lambda_{k+1},\mu_{k+1})
-\nabla_xL(x_k,\lambda_{k+1},\mu_{k+1}).
$$

固定乘子是为了主要捕捉变量移动带来的曲率。BFGS 为

$$
B_{k+1}=B_k-\frac{B_ks_ks_k^TB_k}{s_k^TB_ks_k}
+\frac{y_ky_k^T}{s_k^Ty_k}.
$$

当 $B_k\succ0$ 且 $s_k^Ty_k>0$ 时保持正定。曲率条件失效或分母太小时，应阻尼、跳过或重置；大规模问题可以采用 L-BFGS。正定 BFGS 便于求凸 QP，但不能据此证明真实二阶条件成立。

## 8. 第一篇速查

| 组件 | 回答的问题 |
| --- | --- |
| $J_c$ | 变量移动让约束一阶改变多少？ |
| QP | 如何兼顾局部目标曲率和线性约束？ |
| QP KKT | 如何联立求方向和新的乘子估计？ |
| KKT–Newton | 为什么 QP 是原一阶条件的局部解法？ |
| Globalization | 候选方向可以走多远、能否被接受？ |

**不要误读**：完整步不一定可行；QP 乘子不等于乘子增量；$B_k$ 不只是 $\nabla^2f$；KKT 矩阵不因 $B_k$ 正定就整体正定。

下一篇用可行曲线和三次约束实例解释：为什么真正的二阶结构是 $\nabla^2_{xx}L$，以及如何只在允许的方向上观察它。
