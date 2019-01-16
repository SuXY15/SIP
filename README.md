目录：

[TOC]

### 〇、引言

+ 组织
  - definition of supernova
  - formation of type Ia supernova (some assumption)
  - importance of simulation considering the shape and width of flame area(reacting zone).
+ 方法
  + attemping methods: by fluidic equation considering relativistic effect, considering the influence of large Lewis Number and the acceleration of flamilar and turbulent flame.

1. 关于超新星的演化，前人做出了许多研究工作。Wheeler 和 Harkness 对超新星中的 Ia 型超新星给出了定义：“光谱中没有氢线，并且在光强最大值附近有天鹅 P 型谱线”[^1]
2. I型超新星中的碳氧反应已被化学动力学详细求解过，但基于流体力学的控制方程和状态方程，考虑流场特性及其影响的计算较少。
3. 考虑火焰厚度及大Lewis数下的曲率效应，可以得到较为可观的火焰加速，考虑湍流的强化作用，可以加速到接近声速。

### 一、基本方程

+ 控制方程(可压缩，非定常，无黏)

    $$
    \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho u) = 0
    $$

    $$
    \frac{\partial(\rho u)}{\partial t} + \nabla \cdot (\rho u \otimes u ) + \nabla P= 0
    $$

    $$
    \frac{\partial(\rho E_s)}{\partial t} + \nabla \cdot (\rho u E_s + u P) = \nabla \cdot (\lambda \nabla T)
    $$

    $$
    \frac{\partial(\rho Y)}{\partial t} + \nabla \cdot (\rho u Y) = \nabla \cdot (\rho D \nabla Y) - \omega
    $$

+ 而球坐标下满足：[^1]
    $$
    \nabla F = \frac{\partial F}{\partial r} \vec e_r + \frac{1}{r} \frac{\partial F}{\partial \theta} \vec e_{\theta} + \frac{1}{r\sin \theta}\frac{\partial F}{\partial \varphi} \vec e_{\varphi}
    $$

    $$
    \nabla \cdot F = \frac{1}{r^2}\frac{\partial(r^2 F_r)}{\partial r} + \frac{1}{r\sin\theta} \frac{\partial(F_{\theta} \sin\theta)}{\partial \theta} + \frac{1}{r\sin \theta}\frac{\partial F_{\varphi}}{\partial \varphi}​
    $$

+ 则有一维、球形坐标下的控制方程为：
    * Continuity Equation:
      $$
      \frac{\partial \rho}{\partial t}+\frac{\partial (\rho u)}{\partial r}+\frac{2(\rho u)}{r}=0
      $$

    * Momentum Equation:
      $$
      \frac{\partial (\rho u)}{\partial t}+\frac{\partial (\rho u^2+P)}{\partial r}+\frac{2(\rho u^2)}{r}=0
      $$

    * Energy Equation:
      $$
      \frac{\partial (\rho E_s)}{\partial t}+\frac{\partial (\rho u E_s +u P)}{\partial r}+\frac{2(\rho u E_s + u P)}{r}= \lambda\frac{\partial^2T}{\partial r^2}+\frac{2\lambda}{r}\frac{\partial T}{\partial r}
      $$

    * Convection Diffusion Equation:
      $$
      \frac{\partial (\rho Y)}{\partial t}+\frac{\partial (\rho u Y)}{\partial r} + \frac{2(\rho u Y)}{r}=\frac{\partial}{\partial r}(\rho D \frac{\partial Y}{\partial r})+\frac{2\rho D}{r}\frac{\partial Y}{\partial r}-\omega
      $$

    * Equation of Reaction:
      $$
      \omega = A \rho^k Y^2 exp({\frac{-E_a}{T_9^{\frac{1}{3}}}})
      $$

    * Equation of State:

      $$
      E_{(\rho,T)}=\frac{3}{4{\rho}}(3\pi2){1/3} \hbar c_1({\rho}N)^{4/3} + \frac{1}{2} N \frac{(3\pi2){2/3}}{3 \hbar c_1}(\frac{1}{{\rho} N})^{1/3} (k{T})^2
      $$
      $$ E_s = E_{(\rho,T)}+\frac{1}{2}u^2+qY$$

      $$ P = \frac{\rho E_{(\rho,T)}}{3} $$

      记 $$\gamma = \frac{4}{3} $$，则$$ P=(\gamma -1) \rho E_{(\rho,T)} $$，$$ H=\gamma \rho E_{(\rho,T)} $$

    * 相关数据：[^2][^3][^4]

      反应指前因子：$$ A = 4.26 \times 10^{26} $$ (not sure) 

      反应活化能：$$E_a = 84.165$$

      反应热值：$$ q  = 5.6×10^{17} erg/g = 5.6 × 10^{13} J/kg $$

      密度：$$\rho = 3.5 \times 10^{10} kg / m^3 $$

      层流火焰速度：$$ S_L = 466 m/s $$

      火焰厚度：$$ \delta_f = 0.9 mm = 9 \times 10^{-4} m $$

    * 一些关系：
      $$ \delta_f = \alpha/S_L = \lambda/(\rho_0 C_p S_L) $$

      $$ \omega_0 = \rho_0 S_L Y_0 / \delta_f $$

      $$ Le = \alpha / D $$

+ 解特征声速：

    + 取:$$ U = [\rho, \rho u, \rho E_s, \rho Y]  = [m_1, m_2, m_3, m_4] $$

    + 则 $$ F = [\rho u, \rho u^2 + P, \rho u E_s + u P, \rho u Y] = [m_2, m_2^2/m_1+P，(m_3 + P) m_2 / m_1, m_4 m_2/m_1 ] $$

    + 其中 $$ P = \frac{1}{3}(\rho E_s - \frac{1}{2}\rho u^2 - \rho q Y) = \frac{1}{3}(m3-m_2^2/(2 m_1)-q m4) $$

    + 求解 Jaccobi 矩阵并求得特征值为：
      $$
      \Lambda = [u,  u，u+c, u-c]，c=\sqrt{\gamma P/\rho}
      $$




### 二、无量纲化

+ 以初始状态的特征声速为基准: (\tilde)
    无量纲化                                                        |  基准
     :-------------------------------------------------------------| :---------------------
    $ \rho = \frac{\hat{\rho}}{\hat{\rho_0}} $                     | $ \hat{\rho_0} = 3.5*10^{10} kg/m^3 $ 
    $ u = \frac{\hat u}{\hat{u_0}} $                               | $ \hat{u_0} = c_0 = \sqrt{\gamma P_0 / \rho_0}= 4.6399 \times 10^6 m/s $ 
    $ P = \frac{\hat P}{\hat{\rho_0} \hat{u_0}^2}$                 | $ \hat{\rho_0} \hat{u_0}^2 = 5.6513 \times 10^{23} kg/(m \cdot s^2) [J/m^3] $
    $ E = \frac{\hat E_{\rho,T}}{\hat{u_0}^2} $                    | $ \hat{u_0}^2 = \gamma P_0 / \rho_0 = \gamma (\gamma - 1) E_{(\rho_0, T_0)} = 2.1529 \times 10^{13}  m^2/s^2 [J/kg]$ 
    $ Y = \frac{\hat Y}{\hat{Y_0}} $                               | $ \hat{Y_0} = 1 $ 
    $ q = \frac{\hat q}{\hat q_0} $                                | $ \hat q_0 = \hat{u_0^2} = 2.1529 \times 10^{13}  m^2/s^2 [J/kg] $
    $ r = \frac{\hat r}{\hat {r_0}} $                              | $ \hat {r_0} = \delta_f = 9*10^{-4} m $
    $ t = \frac{\hat t}{\hat{t_0}} $                               | $ \hat{t_0} =\hat {r_0}/ \hat {u_0} = 1.9397 \times 10^{-10} s $ 
    $ \omega = \frac{\hat \omega}{\hat {\omega_0}} $               | $ \hat {\omega_0} = \frac{\hat{\rho_0}\hat{Y_0}S_L}{\delta_f}= 1.8122 \times 10^{16} kg/(m^3 \cdot s)$ 

+ 考虑 $$\omega_0 =\frac{\hat{\rho_0}\hat{Y_0}S_L}{\delta_f} = A \rho_0^k Y_0^2 \exp(-E_a/(T_b/10^9)^{1/3}) $$ ，得指前因子 $$A = \frac{S_L}{\delta_f \hat \rho_0^{k-1} \hat Y_0} \exp(E_a/(T_b/10^9)^{1/3}) = 2.9521 \times 10^{30} \rho_0^{1-k}$$，当前计算取$k=1$

+ 假设反应前后的温度变化与对应焓的变化相同，则：$$c_p(T_b-T_0) = \gamma (E_{(\rho_b,T_b)}-E_{(\rho_0,T_0)})$$，取前后密度相等 $$\rho_b=\rho_0$$，可得到 $$c_p=\frac{\gamma}{\gamma-1} (P_b-P_0)/[\rho_0 (T_b-T_0)] = 5811.8 J/(kg \cdot K)  $$

### 三、无量纲化后得到的方程

+ 控制方程
    * 连续方程：
    $$
    \frac{\partial \rho}{\partial t}+\frac{\partial }{\partial r}(\rho u) + \frac{2}{r}(\rho u)=0
    $$
    * 动量方程：
    $$
    \frac{\partial (\rho u)}{\partial t}+\frac{\partial}{\partial r}(\rho u^2+P)+\frac{2\rho u^2}{r}= 0
    $$
    * 能量方程：
    $$
    \frac{\partial (\rho E_s)}{\partial t}+\frac{\partial}{\partial r}(\rho u E_s + uP) + \frac{2}{r}(\rho u E_s + uP)= A_1(\frac{\partial^2T}{\partial r^2}+\frac{2}{r}\frac{\partial T}{\partial r})
    $$
    * 组分方程：
    $$
    \frac{\partial (\rho Y)}{\partial t}+\frac{\partial}{\partial r}(\rho u Y)+\frac{2}{r}(\rho u Y)=\frac{1}{Le_0C_t} \frac{\partial}{\partial r}(\rho \frac{\partial Y}{\partial r})+\frac{2\rho}{r Le_0C_t}\frac{\partial Y}{\partial r}-\frac{1}{C_t}\omega
    $$
    * 反应方程：
    $$
    \omega = \rho Y \exp(\frac{E_a}{{{T_{b9}}^{1/3}}} - {\frac{E_a}{T_9^{1/3}}} )
    $$
    * 状态方程(能量):
      $$ P = (\gamma -1)\rho E_{(\rho,T)} = (\gamma-1)(\rho E_s - \frac{1}{2} \rho u^2-qcon \rho Y)$$
      $$ E_s =P/(\gamma-1)+\frac{1}{2}u^2+qconY $$

    * 参数：

      $C_t=\frac{\tilde{u_0}}{S_L}$

      $ A_1 = \frac{\lambda(\hat {T_b}-\hat{T_0})}{\hat \rho_0 \hat u_0^3 \hat r_0} = \frac{c_p (T_b-T_0)}{\hat u_0^2 C_t} = 8.4048 \times 10^{-5}$

+ 程序计算
    $$
    U = (\rho,\rho u,\rho E_s,\rho Y)
    $$
    $$
    F = (\rho u, \rho u^2 + P,\rho u E_s + uP,\rho u Y)
    $$
    $$
    G = (\rho u, \rho u^2,\rho u E_s + uP,\rho u Y)
    $$
    $$
    D = (0,0,A_1\frac{\partial^2T}{\partial r^2},\frac{1}{Le_0 C_t} \frac{\partial}{\partial r}(\rho \frac{\partial Y}{\partial r}))
    $$
    $$
    E = (0,0,A_1\frac{\partial T}{\partial r},\frac{\rho}{Le_0 C_t}\frac{\partial Y}{\partial r})
    $$
    $$
    S = (0,0,0,-\frac{\omega}{C_t})
    $$

    - 不考虑曲率：
      $$
      \frac{\partial U}{\partial t}+\frac{\partial F}{\partial r} = D + S
      $$

    - 考虑曲率：
      $$
      \frac{\partial U}{\partial t}+\frac{\partial F}{\partial r} + \frac{2G}{r} = D + \frac{2E}{r} +S
      $$





### 四、问题汇总

#### 1. 偶尔出现数值波动：

- 可能的原因：浮点数计算精度有限

#### 2. 边界条件不对

+ 对称

+ 壁面


#### 3. `ROE`速度为0的处理

+ 特征值中含有速度分母项，无法设置速度为0，否则出现问题

#### 3. 无法加曲率

- 加上曲率后内部速度增加，不满足物理关系
- 对动量方程，$\frac{\partial (\rho u)}{\partial t}+\frac{\partial (\rho u^2+P)}{\partial r}+\frac{2(\rho u^2)}{r}=0$ ，只要速度不为零，必然会不断加速，且越靠近球心，加速越快。

### 五、参考资料

+ [Wikiwand:欧拉方程_(流体动力学)](http://www.wikiwand.com/zh-cn/%E6%AC%A7%E6%8B%89%E6%96%B9%E7%A8%8B_(%E6%B5%81%E4%BD%93%E5%8A%A8%E5%8A%9B%E5%AD%A6))

[^1]: Wheeler J C, Harkness R P, Rep. Prog. Phys. 1990, 53:1467-1557

[^2]: 维基百科球坐标散度 https://en.wikipedia.org/wiki/Divergence#Spherical_coordinates

[^3]: Woosely. 2011. FLAMES IN TYPE Ia SUPERNOVA: DEFLAGRATION–DETONATION TRANSITION IN THE OXYGEN-BURNING FLAME

[^4]: Fowler, W. A., Caughlan, G. R., & Zimmerman, B. A. 1975, ARA&A, 13, 69