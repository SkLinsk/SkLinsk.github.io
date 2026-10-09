# Funnel 与完整 SQP：漏斗边界、恢复机制与复习速查

> **核心结论**：Funnel 维护允许的最大违反度 $\tau_k$，并逐步收紧它。在边界内仍须证明目标或可行性进展；$\tau_k$ 单调不增，不代表 $h_k$ 每步单调下降。

本篇承接 Filter 的步类型与 switching，再给出可逐行复习的 SQP 教学框架。不同文献的接受条件并不完全相同，下面会固定预测量和更新约定，避免混用。

## 1. 从历史边界到漏斗上界

Filter 保存多个 $(h_i,f_i)$；Funnel 主要维护一个违反度上界：

$$
h(x_k)\le\tau_k,\qquad \tau_{k+1}\le\tau_k.
$$

在 $(h,f)$ 平面上，$h>\tau_k$ 的点位于允许区之外。位于边界内只是有资格接受，还要检查相应充分下降条件。

漏斗宽度单调本身不保证 $h\to0$：它可能收敛到正数。完整收敛分析还依赖 switching、目标有界、恢复机制与子问题性质等条件。[统一 funnel restoration SQP 框架](https://arxiv.org/abs/2409.09208)

## 2. 固定一个线搜索教学版本

为与 Armijo 一致，定义 $D_k=-q_k^Td_k$、$P_k(\alpha)=\alpha D_k$。先要求 $h_t\le\tau_k$。

若

$$
D_k>0,\qquad P_k(\alpha)\ge\delta h_k^2,
$$

则按 f-type 验收：

$$
f_t\le f_k-\sigma P_k(\alpha),\qquad \tau_{k+1}=\tau_k.
$$

若不触发该分支，按 h-type 要求

$$
h_t\le\beta\tau_k,\qquad 0<\beta<1,
$$

接受后更新

$$
\tau_{k+1}=(1-\kappa)h_t+\kappa\tau_k,\qquad 0<\kappa<1.
$$

当 $\tau_k>0$ 时，有

$$
h_t\le\tau_{k+1}\le[\kappa+(1-\kappa)\beta]\tau_k<\tau_k.
$$

这既让新点留在更新后的漏斗内，又收窄了上界。某些论文使用二次模型下降 $\Delta m_k$、不同的违反度或其他步类型；替换时必须连同 switching、充分下降与收敛假设一起替换。

## 3. 为什么允许违反度暂时上升？

设 $\tau_0=1$、$h_0=0.7$，一个 h-type 候选为 $h_t=0.4$，取 $\beta=0.8$、$\kappa=0.5$，得到 $\tau_1=0.7$。

下一步若 $h:0.4\to0.5$，仍满足 $h_t<0.7$；只要 f-type 的预测与实际充分下降都通过，就可接受，宽度保持 $0.7$。

| 状态 | 违反度 | 漏斗宽度 | 解释 |
| --- | --- | --- | --- |
| 初始 | $0.7$ | $1.0$ | 位于初始漏斗内 |
| h-type 后 | $0.4$ | $0.7$ | 收窄允许边界 |
| f-type 后 | $0.5$ | $0.7$ | 为目标改善使用边界内的余量 |

这里“h-type”要求相对宽度足够深入，并不严格要求比当前 $h_k$ 更小。Filter 与 Funnel 都不应只凭名称推出逐步单调性。

## 4. Merit / Filter / Funnel 对比

| 方法 | 主要状态 | 对候选点的要求 | 维护任务 |
| --- | --- | --- | --- |
| Merit | 罚参数 $\rho$ | $f+\rho h$ 充分下降 | 处理尺度与罚参数 |
| Filter | 历史 pair 集合 | 对每个条目改善一维，配合当前点测试和 switching | 插入、删除支配条目 |
| Funnel | 违反度上界 $\tau_k$ | 保持在边界内，再满足目标下降或收窄条件 | 保持与更新上界 |

可记成：Merit 设兑换率，Filter 记历史成绩，Funnel 设可行性边界。在某些框架中可用特殊 filter 解释 funnel，但不是任意 funnel 都等于“只存一个历史点的 filter”。

三者都可配合 SOC 与 restoration；线搜索和信赖域则提供拒绝后的步骤控制。

## 5. 停止条件：要检查完整 NLP KKT

采用 $g\le0$、$\mu\ge0$ 的约定，可分别检查

$$
\begin{aligned}
r_s&=\|\nabla_xL(x,\lambda,\mu)\|_\infty,\\
r_e&=\|c(x)\|_\infty,\qquad r_i=\|[g(x)]_+\|_\infty,\\
r_d&=\|[-\mu]_+\|_\infty,\qquad r_c=\|\mu\odot g(x)\|_\infty.
\end{aligned}
$$

各项达到合理缩放后的容差，才报告近似一阶 KKT 收敛。只有等式时，无需不等式对应项。$L$ 的数值本身不是 stationarity residual。

例如 $f(x)=10^{-12}x^2$ 在 $x=10^5$ 时梯度仅为 $2\times10^{-7}$，可能通过宽松绝对容差；因此小残差必须结合尺度解释，不能只看打印出的科学记数法。

步很小或目标变化很小只能作为停滞诊断，不能替代 KKT。小 KKT 残差也不保证局部极小，更不保证一般非凸问题的全局最优。

## 6. 完整 SQP 教学伪代码

这是**模块化线搜索教学框架**，不是某篇论文或生产求解器的逐行实现。Filter 固定采用第三篇的方向导数 switching；Funnel 固定采用本篇的线性预测 $\alpha D_k$。所有恢复退出、曲率修正和 SOC 大小界都须按所选完整算法具体化。

参数：$0<\sigma,r,\beta,\kappa,\gamma_h<1$；$\delta,\gamma_f>0$；步长下界 $\alpha_{\min}>0$。

```text
输入：x0, λ0, μ0≥0, B0, 缩放与容差, 迭代/时间上限
选择接受策略 mode ∈ {MERIT, FILTER, FUNNEL}
初始化 ρ；F=空集；τ0>0 且 τ0≥h(x0)，保留边界余量

for k = 0, 1, ...:
    计算 f, q=∇f, c, g, A=Jc, C=Jg
    若完整 NLP KKT 残差达标：
        返回“一阶 KKT 收敛”及各项残差
    若超过迭代/时间上限：返回“达到上限”

    exact：B=∇²xx L(xk,λk,μk)
    quasi-Newton：使用当前 Bk
    必要时修正曲率或正则化，适配所用 QP 求解策略

    求 QP：min qᵀd + 0.5 dᵀBd
            s.t. c+Ad=0, g+Cd≤0
    得到 d, λhat, μhat
    若 QP 没有有效解：诊断并尝试修复，失败转 RESTORE
    若 d 很小但 KKT 残差仍大：诊断停滞，必要时 RESTORE

    MERIT：调整 ρ，使 φρ'(xk;d)<0
           若无法取得下降方向：RESTORE
    D=-qᵀd；α=1；accepted=false

    while α≥αmin:
        xt=xk+αd，计算真实 ft, ht
        (ok,type)=ACCEPT(xt,α,D,mode)
        若 ok：accepted=true；break

        若 α=1 且拒绝主要源于二阶约束偏差：
            尝试有限次 SOC，求足够小的 w
            xsoc=xk+d+w
            若 xsoc 通过同一真实接受规则：
                xt=xsoc；记录 SOC；accepted=true；break
        α=rα

    若未 accepted：转 RESTORE

    FILTER 且 type=h：将旧点 (hk,fk) 插入 F，删除被支配条目
    FUNNEL 且 type=h：τnew=(1-κ)ht+κτk
    FUNNEL 且 type=f：τnew=τk

    xnew=xt
    选 αdual∈(0,1]，更新乘子估计：
        λnew=λk+αdual(λhat-λk)
        μnew=μk+αdual(μhat-μk)
    或在新点重新估计乘子，保持 μnew≥0

    quasi-Newton：
        s=xnew-xk，包含 SOC 的实际位移
        y=∇xL(xnew,λnew,μnew)-∇xL(xk,λnew,μnew)
        检查曲率与分母，做阻尼 BFGS、跳过或重置
    更新状态并继续

RESTORE:
    解局部违反度最小化或 elastic 子问题
    寻找有足够进展且满足主算法重返条件的 xR
    若成功：
        按算法规则维护 F 或 τ，保证相应不变量
        更新点、重新估计乘子、修复/重置 B，返回主迭代
    否则：报告“恢复失败/不可行驻点/数值失败”等诊断
```

### 6.1 接受模块

```text
ACCEPT(xt,α,D,mode):
    MERIT:
        返回 φρ(xt)≤φρ(xk)+σ α φρ'(xk;d)

    FILTER:
        若 ht 超过允许上限：拒绝
        对 F 中每个 (hi,fi)：
            若 ht≤(1-γh)hi 与 ft≤fi-γf hi 都不成立：拒绝
        switch=(hk≤hmin 且 D>0 且 α D^sf>δ hk^sh)
        若 switch：返回 (ft≤fk-σ αD, f)
        否则返回 (ht≤(1-γh)hk 或 ft≤fk-γf hk, h)

    FUNNEL:
        若 ht>τk：拒绝
        switch=(D>0 且 αD≥δ hk²)
        若 switch：返回 (ft≤fk-σ αD, f)
        否则返回 (ht≤βτk, h)
```

SOC 使用基步的预测量并限制修正大小；具体可证明算法还须明确专用验收细节。Restoration 不能靠随意清空 filter 或放宽 $\tau$ 来绕过进展要求。若 QP 非凸，不能只接受任意驻点方向。

## 7. 复习时最容易混淆的组件

| 组件 | 职责 | 不保证什么 |
| --- | --- | --- |
| SQP / QP | 产生局部候选方向 | 真实非线性可行、完整步可接受 |
| Switching | 选择目标验收或其他进展分支 | 逃离局部极小、全球搜索 |
| SOC | 修补小的二阶约束偏差 | 自动可行或自动通过验收 |
| Restoration | 正常步骤失效时优先降低违反度 | 失败即证明全局无可行解 |
| Reduced Hessian | 检查允许方向的曲率 | 近似矩阵正定就是真实二阶证明 |

## 8. 一页式公式速查

### 模型与方向

$$
L=f+\lambda^Tc+\mu^Tg,\quad \mu\ge0,\quad B\approx\nabla^2_{xx}L.
$$

$$
\min_d q^Td+\tfrac12d^TBd,\qquad c+Ad=0,\quad g+Cd\le0.
$$

### 几何与曲率

$$
Ad_N=-c,\quad AZ=0,\quad d=d_N+Zp,\quad H_R=Z^TBZ.
$$

等式 KKT 点的真实 $H_R\succ0$ 支持严格局部极小；不等式检查临界锥。BFGS 梯度差在新旧两点使用同一组乘子。

### 三种验收与两种补救

**Merit**：$f+\rho h$ 充分下降。**Filter**：对每个历史 pair 至少改善一维，配合 switching / Armijo。**Funnel**：$h_t\le\tau$，f-type 验目标，h-type 收紧上界。

**SOC**：$Aw\approx-c(x+d)$，用小修正消除二阶偏差。**Restoration**：正常优化失效后优先恢复可行性，再满足主算法重返条件。

### 总流程与停止含义

**求导 → 构造 $B$ → 解 QP → 验收 / 缩步 → SOC 或恢复 → 更新点和乘子 → BFGS → 检查完整 KKT。**

小残差是近似一阶条件，真实切向曲率用于局部性质判断；非凸 NLP 一般不保证全局最优。$f$、$h$ 不必每步单调，$\tau$ 单调本身也不保证 $h\to0$。

## 参考与阅读方式

本系列根据提供的 SQP 学习总结重组，保留单位圆与三次约束推导，统一乘子符号，并补充适用条件、数值复核和接受流程。各知识点按“概念 → 数学结构 → 容易卡住的问题”组织。

- [Filter line-search 的局部收敛分析](https://doi.org/10.1137/S1052623403426544)：进一步理解 SOC 与局部收敛。
- [Wächter–Biegler 的实现说明](https://doi.org/10.1007/s10107-004-0559-y)：进一步理解 filter、switching 和 restoration；其算法本体为内点法。
- [A Unified Funnel Restoration SQP Algorithm](https://arxiv.org/abs/2409.09208)：进一步理解漏斗与恢复框架。

不同论文的违反度、预测下降、switching 和历史更新必须一起阅读。本系列伪代码明确标为教学版本，不能直接当作某篇论文的收敛证明对象。
