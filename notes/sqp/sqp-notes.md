# SQP 约束优化：从局部模型到 Filter / Funnel

> **核心结论**：SQP 用一系列二次规划（QP）产生约束优化方向；merit、filter、funnel 决定候选步如何获得接受；SOC 修补非线性约束的二阶偏差；restoration 在正常步骤失效时优先恢复可行性。

**阅读说明**：本文按概念关系重组教学讨论，统一采用“约束 Jacobian 按行排列”的记号。可读取的原对话包含完整 SQP 流程、filter 的步类型与 switching、KKT 残差及 funnel；更早内容未由对话接口返回。因此 Jacobian、三次约束、Hessian、null-space 等部分按本次指定主题补充推导，三次约束例子是补充教学例，非原对话逐字复原。接受条件和伪代码标明其算法版本，不能把不同论文的条件随意拼接成一个有收敛证明的实现。

## 目录

1. 问题、记号与 SQP 的动机
2. Jacobian 与线性化约束
3. QP 子问题、KKT 与 Newton 本质
4. 为什么使用 Lagrangian Hessian
5. Exact SQP 与 quasi-Newton / BFGS
6. 三次约束：一次完整的局部模型计算
7. Reduced Hessian、null-space 与二阶条件
8. Normal / tangential step
9. Globalization 与 merit function
10. Filter-SQP：接受区域与历史更新
11. f-type、h-type 与 switching condition
12. SOC 与 Maratos effect
13. Restoration phase
14. Funnel 方法
15. 完整教学伪代码
16. 常见误区与诊断
17. 网页发布与笔记维护
18. 一页式复习摘要

## 1. 问题、记号与 SQP 的动机

### 1.1 一般非线性规划

考虑光滑问题：

$$
\min_{x\in\mathbb R^n} f(x),\qquad c(x)=0,\qquad g(x)\le0.
$$

其中等式约束有 $m$ 个，不等式约束有 $p$ 个。固定符号约定：

$$
L(x,\lambda,\mu)=f(x)+\lambda^Tc(x)+\mu^Tg(x),\qquad \mu\ge0.
$$

| 符号 | 含义与维数 |
| --- | --- |
| $q_k=\nabla f(x_k)$ | 目标梯度，$n\times1$；避免与不等式函数 $g$ 混淆 |
| $c_k,c(x_k)$ | 等式残差，$m\times1$ |
| $g_k,g(x_k)$ | 不等式函数值，$p\times1$ |
| $A_k=J_c(x_k)$ | 等式 Jacobian，$m\times n$ |
| $C_k=J_g(x_k)$ | 不等式 Jacobian，$p\times n$ |
| $B_k$ | Lagrangian Hessian 或其近似，$n\times n$ |
| $d_k,\alpha_k$ | SQP 方向与步长，$x_{k+1}=x_k+\alpha_kd_k$ |
| $h(x)$ | 约束违反度，非负标量 |

本文一般问题取

$$
h(x)=\|c(x)\|_1+\|[g(x)]_+\|_1,\qquad [g]_j^+=\max(g_j,0).
$$

等式几何推导有时使用二范数；每次涉及数值阈值时应固定范数与缩放。$h$ 的不同定义会改变参数含义。

### 1.2 为什么不是直接沿负梯度走？

负梯度只知道目标在哪里下降，不知道哪些方向会破坏约束。可行域又通常是弯曲的：当前切平面上的直线移动，并不一定停留在真实可行域上。

SQP 把两类信息放进一个局部模型：**用线性约束描述当前允许的移动，用二次模型描述约束优化的曲率**。每次解一个 QP，再在新点重建模型。

> **核心结论**：SQP 的“sequential”指连续重建和求解 QP；不是一次 QP 就能准确替代原非线性问题。

## 2. Jacobian 与线性化约束

### 2.1 $J_c$ 究竟是什么？

将每个约束的梯度转置后堆成行：

$$
J_c(x)=\begin{bmatrix}\nabla c_1(x)^T\\\vdots\\\nabla c_m(x)^T\end{bmatrix}.
$$

Taylor 展开给出

$$
c(x_k+d)=c_k+A_kd+O(\|d\|^2).
$$

因此 QP 要求 $c_k+A_kd=0$，是在要求**线性模型预测下一点可行**。

### 2.2 $A_kd=0$ 什么时候表示切空间？

若 $x_k$ 已经可行且 $A_k$ 满行秩，则可行流形的切空间是

$$
\mathcal T_k=\{d:A_kd=0\}=\ker A_k.
$$

若 $c_k\ne0$，则恢复线性可行性的条件是 $A_kd=-c_k$。这是仿射集合，不是经过原点的切空间。

> **核心结论**：$J_c$ 把“变量怎么动”映射成“约束一阶怎么变”。约束不是一个标量时，不能只拿某一条约束的梯度代替整个 Jacobian。

## 3. QP 子问题、KKT 与 Newton 本质

### 3.1 标准 SQP 子问题

$$
\begin{aligned}
\min_d\quad &q_k^Td+\tfrac12d^TB_kd,\\
\text{s.t.}\quad &c_k+A_kd=0,\\
&g_k+C_kd\le0.
\end{aligned}
$$

常数 $f_k$ 不影响解，故省略。QP 返回方向 $d_k$ 和乘子 $\widehat\lambda_k,\widehat\mu_k$。这里的 QP 乘子是新的乘子估计，**不是默认等于乘子增量**。

QP 的 KKT 条件为

$$
\begin{aligned}
q_k+B_kd_k+A_k^T\widehat\lambda_k+C_k^T\widehat\mu_k&=0,\\
c_k+A_kd_k&=0,\\
g_k+C_kd_k&\le0,\\
\widehat\mu_k&\ge0,\\
\widehat\mu_{k,j}(g_k+C_kd_k)_j&=0.
\end{aligned}
$$

若 QP 凸，这些条件在适当正则性下刻画最优解；若 QP 非凸，满足 KKT 本身不保证它是最小点。

### 3.2 等式情形的 KKT 线性系统

只含等式时：

$$
\begin{bmatrix}B_k&A_k^T\\A_k&0\end{bmatrix}
\begin{bmatrix}d_k\\\widehat\lambda_k\end{bmatrix}
=-\begin{bmatrix}q_k\\c_k\end{bmatrix}.
$$

注意这个块矩阵通常是不定的，即使 $B_k$ 正定，也不能把整个 KKT 矩阵当成正定矩阵。

### 3.3 为什么说 SQP 是对 KKT 做 Newton？

等式 NLP 的方程是

$$
F(x,\lambda)=\begin{bmatrix}\nabla_xL(x,\lambda)\\c(x)\end{bmatrix}=0.
$$

在当前点线性化：

$$
\begin{bmatrix}\nabla^2_{xx}L_k&A_k^T\\A_k&0\end{bmatrix}
\begin{bmatrix}d_k\\\Delta\lambda_k\end{bmatrix}
=-\begin{bmatrix}\nabla_xL_k\\c_k\end{bmatrix}.
$$

用 $\widehat\lambda_k=\lambda_k+\Delta\lambda_k$，第一行就变成 QP stationarity。因此 exact SQP 的等式步骤与 KKT-Newton 步一致。对不等式，局部固定活跃集后有类似解释，但还要处理活跃集变化和互补条件。

## 4. 为什么使用 Lagrangian Hessian

### 4.1 目标曲率之外，还有约束曲率

$$
\nabla^2_{xx}L
=\nabla^2f+\sum_{i=1}^m\lambda_i\nabla^2c_i
+\sum_{j=1}^p\mu_j\nabla^2g_j.
$$

“对 $x$ 求二阶导”时乘子固定。它不是同时对 $x,\lambda,\mu$ 求 Hessian。$B_k$ 中的约束曲率并不意味着 QP 约束也变成二次约束；QP 约束依然线性。

### 4.2 沿可行曲线推导

只考虑等式，在 KKT 点 $x_*$ 上取一条可行曲线 $x(t)$，满足 $x(0)=x_*$、$x'(0)=v$。因为 $c(x(t))=0$，

$$
A_*v=0,\qquad
A_*x''(0)+\begin{bmatrix}v^T\nabla^2c_1v\\\vdots\\v^T\nabla^2c_mv\end{bmatrix}=0.
$$

目标沿曲线的二阶导包含“路径转弯”这一项：

$$
\frac{d^2}{dt^2}f(x(t))\bigg|_0
=v^T\nabla^2fv+\nabla f^Tx''(0).
$$

由 stationarity $\nabla f=-A_*^T\lambda_*$，代入得

$$
\frac{d^2}{dt^2}f(x(t))\bigg|_0
=v^T\left(\nabla^2f+\sum_i\lambda_{*,i}\nabla^2c_i\right)v
=v^T\nabla^2_{xx}Lv.
$$

> **核心结论**：$\nabla^2f$ 看直线上的目标曲率；$\nabla^2_{xx}L$ 把沿弯曲可行域运动的效应计入。这也是 KKT-Newton 系统自然出现它的原因。

## 5. Exact SQP 与 quasi-Newton / BFGS

### 5.1 Exact SQP

每次计算 $B_k=\nabla^2_{xx}L(x_k,\lambda_k,\mu_k)$。在解附近、LICQ、二阶充分条件及足够光滑等假设下，完整 Newton 步可以获得很快的局部收敛。远离解时，真实 Hessian 可能不定，QP 可能无界；需要修正曲率、信赖域或合适的非凸子问题处理。

### 5.2 Quasi-Newton SQP

不用显式二阶导，而利用梯度变化近似 Lagrangian Hessian。令实际接受的位移为

$$
s_k=x_{k+1}-x_k,
$$

在新旧两点使用**相同的更新后乘子**构造

$$
y_k=\nabla_xL(x_{k+1},\lambda_{k+1},\mu_{k+1})
-\nabla_xL(x_k,\lambda_{k+1},\mu_{k+1}).
$$

这样主要捕捉 $x$ 改变带来的曲率，而不是乘子改变本身。BFGS 更新为

$$
B_{k+1}=B_k-\frac{B_ks_ks_k^TB_k}{s_k^TB_ks_k}
+\frac{y_ky_k^T}{s_k^Ty_k}.
$$

当 $B_k\succ0$ 且 $s_k^Ty_k>0$ 时保持正定。约束问题中这一曲率条件可能失败，应阻尼、跳过更新或重置。

一种 Powell 阻尼取 $0<\eta_B<1$：若 $s^Ty<\eta_Bs^TBs$，令

$$
\theta=\frac{(1-\eta_B)s^TBs}{s^TBs-s^Ty},\qquad
\bar y=\theta y+(1-\theta)Bs,
$$

否则 $\bar y=y$；用 $\bar y$ 更新，并检查分母与数值尺度。

| 方法 | 曲率信息 | 主要注意点 |
| --- | --- | --- |
| Exact | 真实 $\nabla^2_{xx}L$ | 二阶导成本、负曲率、正则化 |
| BFGS | 梯度差、secant 近似 | 固定乘子构造 $y$，阻尼与重置 |
| Reduced BFGS | 只近似切空间曲率 | 基底随迭代改变，需要一致的坐标处理 |

正定 BFGS 便于产生凸 QP，但它不能作为“真实二阶条件成立”的证据。

## 6. 三次约束：一次完整的局部模型计算

**补充教学例**：在三次曲线 $y=x^3$ 上最小化

$$
f(x,y)=\tfrac12[(x-1)^2+(y-1)^2],\qquad c(x,y)=y-x^3=0.
$$

目标 Hessian 是 $I$；约束的信息是

$$
J_c(x,y)=\begin{bmatrix}-3x^2&1\end{bmatrix},\qquad
\nabla^2c(x,y)=\begin{bmatrix}-6x&0\\0&0\end{bmatrix}.
$$

因此

$$
\nabla^2_{xx}L=\begin{bmatrix}1-6\lambda x&0\\0&1\end{bmatrix}.
$$

### 6.1 在一个不可行点构造 QP

取 $(x_k,y_k)=(1,0)$、$\lambda_k=-1/6$。则

$$
c_k=-1,\quad q_k=\begin{bmatrix}0\\-1\end{bmatrix},\quad
A_k=\begin{bmatrix}-3&1\end{bmatrix},\quad B_k=\begin{bmatrix}2&0\\0&1\end{bmatrix}.
$$

设 $d=(u,v)$，QP 为

$$
\min_{u,v}\ -v+u^2+\tfrac12v^2,\qquad -1-3u+v=0.
$$

代入 $v=1+3u$，目标化为 $-1/2+(11/2)u^2$。故 $d_k=(0,1)$，QP 乘子为 $0$；完整步到达 $(1,1)$，这是 $f=0$ 的全局最小点。

这次恰好沿 $y$ 修正，非线性约束也精确满足。一般方向没有这个幸运：

$$
c(x+u,y+v)=c(x,y)-3x^2u+v-3xu^2-u^3.
$$

QP 只消掉前三项，留下 $-3xu^2-u^3$。这直接展示了 SOC 为什么有用。

### 6.2 切方向与约束曲率

在可行点 $(a,a^3)$，可以取 null-space 基底

$$
Z=\begin{bmatrix}1\\3a^2\end{bmatrix},\qquad J_cZ=0.
$$

沿该基底的 reduced Hessian 为

$$
Z^T\nabla^2_{xx}LZ=1-6\lambda a+9a^4.
$$

若只用 $\nabla^2f=I$，会漏掉 $-6\lambda a$。基底不必归一化，但比较特征值大小时要注意基底尺度；正定性在可逆换基下不变。

## 7. Reduced Hessian、null-space 与二阶条件

### 7.1 为什么投影到 null-space？

等式正则点的允许一阶方向满足 $Av=0$。若 $Z\in\mathbb R^{n\times(n-m)}$ 的列张成 $\ker A$，则 $v=Zp$，

$$
v^TBv=p^T(Z^TBZ)p.
$$

$Z^TBZ$ 是 reduced Hessian。等式问题在 LICQ 和 KKT 下，若

$$
Z^T\nabla^2_{xx}L_*Z\succ0,
$$

则满足严格局部极小的二阶充分条件。完整 Hessian 不必在整个空间正定；真正关心的是可行切方向上的曲率。

不等式问题的二阶条件应在 **临界锥（critical cone）** 上检查。只有在合适的活跃集、严格互补等条件下，才可简化为活跃约束 Jacobian 的 null-space 检查。

### 7.2 Null-space 解 QP

先找 $d_N$ 满足 $Ad_N=-c$，再写 $d=d_N+Zp$。于是

$$
(Z^TBZ)p=-Z^T(q+Bd_N).
$$

这把约束 QP 的自由度压缩到 $n-m$ 维。若 $A$ 满行秩且 $Z^TBZ\succ0$，等式 QP 有唯一最小解，KKT 矩阵也非奇异。

## 8. Normal / tangential step

### 8.1 两个职责

写成 $d=d_N+d_T$：normal step 主要降低线性化约束残差，tangential step 主要利用余下自由度改善目标。

等式满行秩时，最小范数 normal step 为

$$
d_N=-A^T(AA^T)^{-1}c.
$$

实际计算用 QR、SVD 或线性求解，避免显式形成逆矩阵。Tangential step 满足 $Ad_T=0$，并求解

$$
\min_{d_T}\ (q+Bd_N)^Td_T+\tfrac12d_T^TBd_T,\qquad Ad_T=0.
$$

### 8.2 远离可行域时

在信赖域方法中，normal step 常改为

$$
\min_{\|d_N\|\le\Delta_N}\ \tfrac12\|c+Ad_N\|_2^2.
$$

这允许只做部分可行性修复，再用切向步完成优化。Tangential step 还须遵守总步长预算。对不等式，需要结合活跃集或 slack 变量；不能把所有不等式都当作必须满足等号的约束。

> **核心结论**：normal / tangential 是一种几何和算法分工，不是所有 SQP 实现都必须显式求两个独立子问题。

## 9. Globalization 与 merit function

### 9.1 为什么不能一直走完整步？

局部模型只在附近准确。远离解时，完整步可能使非线性约束恶化、目标上升或进入模型失真的区域。Globalization 用线搜索或信赖域控制这一风险。

这里“global convergence”通常指从较广的初始点范围建立到驻点等的收敛性质，**不是保证找到全局最优解**。

### 9.2 Merit：把两个目标合成一个标量

常用 $\ell_1$ 精确罚函数：

$$
\phi_\rho(x)=f(x)+\rho h(x),\qquad \rho>0.
$$

沿方向要求 Armijo：

$$
\phi_\rho(x_k+\alpha d_k)
\le\phi_\rho(x_k)+\sigma\alpha\phi_\rho'(x_k;d_k),\qquad 0<\sigma<1.
$$

因为 $h$ 不光滑，这里是单侧方向导数。只有方向导数为负时，该式才要求下降；否则需调整罚参数或修复方向。

等式 QP 满足 $Ad=-c$ 时，$\|c\|_1$ 的方向导数是 $-\|c\|_1$。在适当条件下，罚参数大于解处乘子的相应对偶范数，可体现局部“精确罚”性质；但参数过大容易让优化被可行性项压住。

## 10. Filter-SQP：接受区域与历史更新

### 10.1 Filter 保存什么？

保存若干历史 pair：$\mathcal F=\{(h_i,f_i)\}$。低 $h$、低 $f$ 都更好。候选点必须对**每一个** filter 条目至少在一个维度取得带余量的改善，例如

$$
\forall(h_i,f_i)\in\mathcal F:\quad
h_t\le(1-\gamma_h)h_i\quad\text{或}\quad
f_t\le f_i-\gamma_fh_i.
$$

其中 $0<\gamma_h<1$、$\gamma_f>0$。这里的“或”在每个条目内部成立；不是找到任意一个条目可以改善就通过。

若新 pair 两个坐标都不比某条目更大，则新 pair 支配该条目，后者可删除。加入历史点使算法不容易反复回到已经取得过的较差折中。

### 10.2 只有 filter 检查够吗？

不够。还要与当前点比较、设置违反度上限，并用 switching 和目标充分下降条件防止接近可行域时只靠微小的约束改善停滞。

> **核心结论**：filter 允许目标暂时上升，也允许违反度暂时上升，但两者不是任意上升。它管理的是经过余量保护的历史接受区域。[Filter line-search 的原始分析](https://doi.org/10.1137/S1052623403426544)

## 11. f-type、h-type 与 switching condition

### 11.1 步类型是验收职责

- **f-type**：方向预测了值得追求的目标下降，必须验证真实目标充分下降。
- **h-type**：没有进入上述目标下降分支，靠 filter 或可行性规则取得接受，目标可以上升。

在某些 filter 算法中，h-type 的逻辑是“非 f-type 的 filter 接受步”，不一定严格等价于 $h_{k+1}<h_k$。名称是帮助理解职责，最终以接受条件为准。

例如 $(f,h):(10,0.5)\to(10.4,0.05)$ 可作为可行性进展；若 $h\approx10^{-7}$，只改善一点 $h$ 而完全不改善最优性，就不能长期满足算法的收敛需求。

### 11.2 简化的模型下降 switching

定义 SQP 局部模型下降：

$$
\Delta m_k(d)=-q_k^Td-\tfrac12d^TB_kd.
$$

它是局部 QP 模型的预测量；因 $B_k$ 含约束曲率，不能把它当作真实 $f$ 的精确 Taylor 下降。一种教学规则是

$$
\Delta m_k(d_k)>0,\qquad \Delta m_k(d_k)\ge\delta h_k^2.
$$

$h$ 小时阈值小，算法更容易要求目标下降。指数 $2$ 和参数并非所有 filter / funnel 方法的统一标准。

### 11.3 本笔记伪代码采用的线搜索规则

为让方向导数和 Armijo 一致，定义 $D_k=-q_k^Td_k$。Filter 分支采用一种 Wächter–Biegler 风格的规则：

$$
h_k\le h_{\min},\qquad D_k>0,\qquad
\alpha D_k^{s_f}>\delta h_k^{s_h}.
$$

触发时用

$$
f(x_k+\alpha d_k)\le f_k-\sigma\alpha D_k.
$$

例如 $f_k=10$、$D_k=2$、$\sigma=0.1$、$\alpha=1$，则要求 $f_t\le9.8$；$9.99$ 虽下降，却不够。

未触发时采用与当前点的 filter 式进展检查。具体论文还可能对指数、阈值和 filter 更新施加额外要求。IPOPT 是采用 filter 的**内点算法**，不能把它整体称为经典 SQP；这里借用其接受逻辑帮助阅读。[Wächter–Biegler 算法说明](https://doi.org/10.1007/s10107-004-0559-y)

### 11.4 历史更新的一个明确约定

本笔记在 h-type 接受后，把**接受前当前点** $(h_k,f_k)$ 插入 filter，删除被该 pair 支配的旧条目；f-type 不插入。不同论文可能采用移动后的点或其他更新时机，应连同接受规则整体阅读。

## 12. SOC 与 Maratos effect

### 12.1 线性可行不等于非线性可行

若 $c_k+Ad=0$，仍可能有

$$
c(x_k+d)=O(\|d\|^2).
$$

Maratos effect 指解附近本来适合快速收敛的完整 SQP 步，因为这些二阶约束误差而被 merit / 接受规则拒绝；反复缩步会妨碍预期的快速局部收敛。

### 12.2 Second-order correction

在完整步 $x_t=x_k+d$ 被约束误差阻挡时，求一个小修正 $w$：

$$
A_kw\approx-c(x_t),\qquad x_{\mathrm{soc}}=x_k+d+w.
$$

也可用更新后的 Jacobian；在正则情形下 $w=O(\|d\|^2)$。修正后仍须重新检查真实函数值和约束值。处理不等式时主要修正有关活跃约束，并验证所有不等式。

### 12.3 单位圆上的直观例子

在 $(1,0)$、约束 $x^2+y^2-1=0$ 上，切向步 $d=(0,t)$ 满足 $Ad=0$，却产生残差 $t^2$。取 $w=(-t^2/2,0)$，有

$$
c(1-t^2/2,t)=t^4/4.
$$

一个二阶小修正将违反度从二阶降为四阶。实际最优化中的 Maratos effect 还涉及目标和接受测试；此例只展示约束修正机制。

> **核心结论**：SOC 是对一个已有好方向的小修补，不是另一次完整优化；它不保证一定可接受，也不是每次回溯都必须调用。

## 13. Restoration phase

### 13.1 什么时候进入？

线性化约束不相容、QP 求解失败且修正无效，或回溯到最小步长仍无法接受时，可以进入 restoration。QP 失败也可能来自数值病态、无界或求解器问题，应先辨别原因。

等式违反度子问题常取

$$
\min_x\ v(x)=\tfrac12\|c(x)\|_2^2,
$$

一般约束可加入 $\tfrac12\|[g(x)]_+\|_2^2$，也可用 slack / elastic 变量、邻近项或信赖域。

### 13.2 退出并不是“$h$ 降了就回去”

必须满足主算法的重返条件：例如获得足够可行性进展，且新点能通过 filter 或进入 funnel；还要能够重新构造可用的主问题步骤。Restoration 后乘子应重新估计，Hessian 近似应按规则重置或修复。

若 restoration 停在 $h>0$ 的驻点，只能报告恢复失败或不可行驻点等状态，**不能仅凭它断言原问题无可行解**。

> **核心结论**：restoration 是优先恢复可行性的独立阶段；normal step 是正常步骤内部的组成部分，二者不是同义词。

## 14. Funnel 方法

### 14.1 从二维历史到一个违反度上界

维护 funnel width $\tau_k$，并保证

$$
h(x_k)\le\tau_k,\qquad \tau_{k+1}\le\tau_k.
$$

候选点首先必须处于漏斗内。Funnel 边界只是资格条件，还要检查目标下降或可行性进展。下面给出与原讨论相近的教学规则；完整算法的收敛依赖 switching、目标有界和恢复机制等额外条件，单调 $\tau$ 本身不保证 $h\to0$。[统一 funnel restoration SQP 框架](https://arxiv.org/abs/2409.09208)

### 14.2 f-type 与 h-type

模型版本可比较 $\Delta m_k$ 与 $\delta h_k^2$，要求模型下降为正。若触发 f-type 且

$$
h_t\le\tau_k,\qquad f_k-f_t\ge\sigma\Delta m_k,
$$

则接受并令 $\tau_{k+1}=\tau_k$。有些方法用线性预测下降而非二次模型，必须与实际定义一致。

若未触发，则 h-type 的一种接受条件和更新为

$$
h_t\le\beta\tau_k,\qquad
\tau_{k+1}=(1-\kappa)h_t+\kappa\tau_k,
\qquad 0<\beta,\kappa<1.
$$

于是 $h_t\le\tau_{k+1}<\tau_k$（当 $\tau_k>0$）。它保证收窄上界，但不必保证 $h_t<h_k$。

### 14.3 数值例子

设 $\tau_0=1$、$h_0=0.7$，h-type 候选 $h_t=0.4$，$\beta=0.8$、$\kappa=0.5$。接受后 $\tau_1=0.7$。

下一次 f-type 将 $h:0.4\to0.5$，仍在 $\tau_1=0.7$ 内；若目标充分下降，仍可接受，宽度保持 $0.7$。

> **核心结论**：funnel width 单调不增，约束违反度不必每一步单调不增。某些框架可把 funnel 解释成特殊 filter，但一般 funnel 算法不能不加条件地等同于“只存一个历史点的 filter”。

### 14.4 三种 globalization 的比较

| 方法 | 保存的状态 | 核心接受依据 | 主要维护任务 |
| --- | --- | --- | --- |
| Merit | 罚参数 $\rho$ | $f+\rho h$ 充分下降 | 平衡尺度、更新罚参数 |
| Filter | 若干 $(h_i,f_i)$ | 对历史条目改善至少一维，并配合 switching | 插入、支配删除、历史测试 |
| Funnel | 上界 $\tau_k$ | 留在边界内，配合目标下降或收窄条件 | 维护和收紧违反度上界 |

三者都可以配合 SOC 和 restoration。它们管理迭代的进展，不负责在不同非凸吸引域之间搜索全局最优解。

## 15. 完整教学伪代码

### 15.1 范围与约定

以下是**线搜索 SQP 的模块化教学框架**，覆盖成功、回溯、SOC、恢复与终止路径。它不是某个求解器的逐行复现，也不是直接运行的生产实现。Filter 使用第 11 节的方向导数规则；funnel 为保持线搜索一致，采用 $\alpha D_k$ 的线性预测版本。第 14 节的二次模型版本是另一选择，不与这里混用。

参数满足 $0<\sigma,r,\beta,\kappa,\gamma_h<1$，$\delta,\gamma_f>0$。停止使用经过合理缩放的完整 NLP KKT 残差，包含 stationarity、等式、不等式违反度、乘子非负性和互补性。

```text
输入 x0, λ0, μ0 ≥ 0, B0, 容差, 迭代/时间上限
选择 mode ∈ {MERIT, FILTER, FUNNEL}
初始化罚参数 ρ；filter F = 空集；τ0 ≥ h(x0)

for k = 0, 1, ...:
    计算 f, q=∇f, c, g, A=Jc, C=Jg
    若完整、缩放后的 NLP KKT 残差满足容差：
        返回“一阶 KKT 收敛”及残差（不宣称全局最优）
    若迭代/时间达到上限：返回“达到上限”

    exact: B = ∇²xx L(xk, λk, μk)
    quasi-Newton: 使用已有 Bk
    必要时按所用 QP 策略修正曲率/正则化

    解 QP: min qᵀd + 0.5 dᵀBd
           s.t. c+Ad=0, g+Cd≤0
    得到 d, λhat, μhat
    若 QP 无有效解：转 RESTORE
    若 d 很小但 KKT 残差仍大：诊断/修复；失败则 RESTORE

    若 MERIT：更新 ρ，确保 φρ'(xk;d)<0
               若无法取得下降方向：RESTORE
    D = -qᵀd；α = 1；accepted = false

    while α ≥ αmin:
        xt = xk + αd；计算 ft, ht
        (ok, type) = ACCEPT(xt, α, D, mode)
        若 ok：accepted = true；break

        若 α=1 且拒绝主要来自非线性约束偏差：
            尝试有限次 SOC，求小修正 w
            xsoc = xk + d + w
            若 w 足够小且 xsoc 通过同一真实接受规则：
                xt = xsoc；记录 SOC；accepted = true；break
        α = r α

    若未 accepted：转 RESTORE

    若 FILTER 且 type=h：
        向 F 插入旧点 (hk, fk)，删除被它支配的条目
    若 FUNNEL 且 type=h：τnew=(1-κ)ht+κτk
    若 FUNNEL 且 type=f：τnew=τk

    xnew = xt
    选择乘子更新策略：如用 αdual∈(0,1]
    λnew = λk + αdual(λhat-λk)
    μnew = μk + αdual(μhat-μk)
    或在新点重新估计乘子；保持 μnew≥0
    若 quasi-Newton：
        s = xnew-xk（含 SOC 的实际位移）
        y = ∇xL(xnew,λnew,μnew)-∇xL(xk,λnew,μnew)
        用阻尼 BFGS 更新；异常时跳过或重置
    更新状态；continue

RESTORE:
    解局部违反度最小化/elastic 子问题
    寻找具有足够违反度进展且满足主算法重返条件的 xR
    若成功：
        按该算法规则维护 F 或 τ，使新点满足不变量
        xnew=xR；重新估计乘子；重置/修复 B；continue
    否则：返回“恢复失败/不可行驻点/数值失败”等诊断
```

### 15.2 接受模块

```text
ACCEPT(xt, α, D, mode):
    MERIT:
        返回 φρ(xt) ≤ φρ(xk)+σ α φρ'(xk;d)

    FILTER:
        若 ht 超过允许上限：拒绝
        若对 F 任一条目，两条改善不等式都不满足：拒绝
        switch = (hk≤hmin 且 D>0 且 α D^sf > δ hk^sh)
        若 switch：返回 (ft≤fk-σ αD, f)
        否则：返回
            (ht≤(1-γh)hk 或 ft≤fk-γf hk, h)

    FUNNEL（线性预测的教学版本）:
        若 ht>τk：拒绝
        switch = (D>0 且 αD≥δ hk²)
        若 switch：返回 (ft≤fk-σ αD, f)
        否则：返回 (ht≤βτk, h)
```

SOC 候选采用其基步的预测量并限制修正大小；有收敛证明的实现还需明确 SOC 次数、大小界和专用验收细节。Restoration 的 filter 插入及 funnel 收窄必须与所选择的论文一致，不能用“随意清空历史或放宽边界”替代。

## 16. 常见误区与诊断

### 16.1 KKT 残差很小，目标为何仍然很高？

取单位圆上的问题 $\min x$，约束 $x^2+y^2=1$。

在 $(1,0)$，取 $\lambda=-1/2$，则 $c=0$、$\nabla_xL=0$，但该点是最大点；$Z=(0,1)^T$，有 $Z^T\nabla^2LZ=-1$。

在 $(-1,0)$，取 $\lambda=1/2$，则 reduced Hessian 为 $1$，目标才取得最小值 $-1$。

因此残差小首先意味着近似一阶驻定和可行。病态或退化时，残差小甚至不一定意味着变量在距离上接近某个精确 KKT 点。

### 16.2 复习时最容易混淆的判断

| 误区 | 正确理解 |
| --- | --- |
| $L$ residual 就是 $|L|$ | stationarity residual 是 $\|\nabla_xL\|$ |
| 只检查 $c$ 和 $\nabla L$ 足够 | 有不等式还要检查违反度、对偶可行与互补性 |
| $\nabla f=0$ 才能最优 | 约束最优时梯度可由约束梯度与乘子抵消 |
| $B$ 一定是 $\nabla^2f$ | 它近似 $\nabla^2_{xx}L$，含约束曲率 |
| $B\succ0$ 才能是约束极小点 | 等式问题主要看切空间曲率；不等式看临界锥 |
| QP 可行则新点一定可行 | 只保证线性化约束；仍有高阶误差 |
| h-type 每步必须降低 $h$ | 依具体 filter / funnel 接受规则，不能只凭名称判断 |
| 在 funnel 内就接受 | 还要满足相应充分下降条件 |
| $\tau$ 单调即可推出 $h\to0$ | 还依赖 switching、界限和恢复等机制 |
| Switching 能逃离局部最优 | 它管理步的验收，不是全局搜索 |
| Restoration 失败证明问题不可行 | 可能停在局部违反度驻点或数值失败 |

### 16.3 一个实用的诊断次序

先检查导数和尺度，再检查完整 KKT 残差；若残差小但解不理想，检查真实 reduced Hessian / 临界锥曲率，并区分局部解和全局解。若残差不小但步长不断缩小，查看 QP 状态、预测下降、实际下降、filter / funnel 拒绝原因以及 SOC、restoration 是否正常工作。

例如 $f(x)=10^{-12}x^2$ 在 $x=10^5$ 时梯度仅为 $2\times10^{-7}$，绝对容差 $10^{-6}$ 可能误判 stationarity。应结合缩放、相对尺度和问题物理含义解释数值。

## 17. 网页发布与笔记维护

把 Markdown 作为可编辑主版本，PDF 作为打印和固定版本。此文件使用 `$...$` 行内数学与 `$$...$$` 独立数学；网页需启用兼容的 KaTeX 或 MathJax 渲染。普通 Markdown 渲染器不一定支持数学。

建议每个知识点保持三层：**它解决什么问题 → 关键公式与算法位置 → 自己曾卡住的“为什么”**。例如“为什么用 $\nabla^2L$”放在 Hessian 推导之后，比按聊天顺序保存零散问答更容易复习。

发布前检查公式、矩阵换行、表格在手机上的宽度和代码块换行。保留术语、符号约定和版本说明；遇到新论文时，把该论文的 $h$ 定义、预测下降和 switching 单独记录，避免偷换符号。

## 18. 一页式复习摘要

### 问题与局部方向

$$
\min f(x),\ c=0,\ g\le0;\qquad L=f+\lambda^Tc+\mu^Tg,\ \mu\ge0.
$$

$$
\min_d q^Td+\tfrac12d^TBd,\qquad c+Ad=0,\quad g+Cd\le0.
$$

$A=J_c$ 按行堆约束梯度；$B\approx\nabla^2_{xx}L$ 包含约束曲率。Exact 使用真实二阶导；BFGS 用同一乘子下的 Lagrangian 梯度差并检查曲率。

### 几何与二阶信息

$Ad_N=-c$ 修复线性可行性；$d_T=Zp$、$AZ=0$ 使用切自由度。Reduced Hessian 为 $Z^TBZ$；真实二阶充分条件在等式切空间、不等式临界锥上判断。

### 三种验收思路

**Merit**：$f+\rho h$ 充分下降。**Filter**：对每个历史 pair 至少改善一维，近可行时加 switching / Armijo。**Funnel**：$h_t\le\tau$，f-type 要目标下降，h-type 收窄上界。

### f / h 与两个补救模块

f-type 用正的预测目标下降触发并验真实下降；h-type 按相应可行性 / filter 规则验收。Switching 的指数和预测量依论文变化。SOC 修二阶约束误差；restoration 在正常步骤失效后优先降低违反度。

### 算法骨架与停止含义

**求导 → 构造 / 更新 $B$ → 解 QP → 验收 / 回溯 → SOC 或 restoration → 更新点与乘子 → BFGS → 检查完整 KKT。**

KKT 包括 stationarity、原始可行、对偶可行、互补性。小残差说明近似一阶条件；真实正的切向曲率支持局部极小；非凸问题仍不保证全局最优。$f$ 和 $h$ 都不必逐步单调。

## 参考与进一步阅读

- 原教学对话：《不会嫌弃你》，本笔记按指定主题重组并补充；不保留聊天中的临时引用占位符。
- [Wächter–Biegler：Filter line-search 的局部收敛分析](https://doi.org/10.1137/S1052623403426544)：用于核对 filter 与 SOC 的关系。
- [Wächter–Biegler：内点 filter line-search 实现](https://doi.org/10.1007/s10107-004-0559-y)：用于理解 switching、restoration、SOC；算法本体为内点法。
- [Kiessling、Leyffer、Vanaret：A Unified Funnel Restoration SQP Algorithm](https://arxiv.org/abs/2409.09208)：用于进一步核对 funnel 边界与恢复框架。本文的具体接受模块仍标为教学版本。
