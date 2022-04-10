<h2 align=center>Spherically Inward Propagating Flame of Type Ia Supernova</h2>

[TOC]

### 1. Introduction

Supernova represent the catastrophic explosions that mark the end of the life of some stars. The ejected mass is of order 1 to 10 solar masses with bulk velocities ranging from a few thousand to a few tens of thousands of km/s.[^1] The traditional single Chandrasekhar mass C–O White Dwarf burning is still considered to be responsible for a large population of type Ia supernova (SN Ia). Specifically, one of the key issues in its modeling is related to the flame acceleration and deflagration-detonation transition (DDT), with flame front instabilities being considered as a possible mechanism in driving the acceleration. The perspective of this proposal is to give a set of solutions of this problem. The milestones of this project are planed to be:

​	(a) Complete a set of code that can solve the 1D $H_2$-$O_2$ flame.

​	(b) Complete a set of code that can solve the simplest reacting flow in SN Ia conditions.

​	(c) Try to observe spherical flame acceleration with large Lewis number curvature effect.

​	(d) Try to observe pulsation in Xing & Zhao's case.

​	(e) Try to observe deflagration-detonation transition (DDT) in (d).



### 2. Basic Equations

#### 2.1 Planar flame

+ Governing equations for 1-dimensional, Cartesian coordinates, compressible and non-viscous reacting flow are:
    $$
    \begin{align*}
    \frac{\partial \boldsymbol U}{\partial t} + \frac{\partial \boldsymbol F(\boldsymbol U)}{\partial x} = \frac{\partial \boldsymbol D(\boldsymbol U)}{\partial x} + \boldsymbol S
    \end{align*}
    $$
    where $\boldsymbol U$ is the conserved variables vector, $\boldsymbol F$ is  the flux vector, $\boldsymbol D$ is the diffusion flux, and $\boldsymbol S$ is the source term vector. The vectors are given by:
    $$
    \begin{align*}
    \boldsymbol U = \left[\begin{matrix} \rho \\ \rho u \\ \rho E \\ \rho Y\end{matrix}\right],
    \,\,\,\,
    \boldsymbol F(\boldsymbol U) = \left[\begin{matrix} \rho u \\ \rho u^2+p \\ \rho u(E+p/\rho) \\ \rho uY\end{matrix}\right],
    \,\,\,\,
    \boldsymbol D(\boldsymbol U) = \left[\begin{matrix} 0 \\ 0 \\  \frac{\mu C_p}{Pr} \frac{\partial T}{\partial x} + \frac{\mu Q}{Sc} \frac{\partial Y}{\partial x} \\ \frac{\mu}{Sc} \frac{\partial Y}{\partial x} \end{matrix}\right],
    \,\,\,\,
    \boldsymbol S =\left[\begin{matrix} 0 \\ 0 \\ 0 \\  \rho \omega \end{matrix}\right]
    \end{align*}
    $$
    Here $\rho$ is density, $u$ is velocity, $E$ is the total energy per unit mass and $Y$ is the reactant mass fraction, respectively. The total energy is consisted by internal energy, kinetic energy and the "chemical" energy in reactant:
    $$
    \rho E = \rho e + \frac{1}{2}\rho u^2 \color{red}{+ \rho QY}
    $$
    
    > The viscosity $\mu = \rho \nu$
    >
    > The Schmidt number $Sc=\nu/D = \mu/\rho D $
    >
    > The Prandtl number $Pr=\nu/\alpha= \mu c_p/\lambda$
    >
    > The enthalpy per unit mass $H=E+p/\rho$
    >
    > <p style="color:red">The heat value part might need to be put in internal energy part</color>
    
    + **In the scenario of methane flame**, the internal energy per unit mass is given by:
      $$
      e = \frac{p}{\rho(\gamma-1)} = \frac{R_u T}{\gamma-1}
      $$
      and the reaction rate is
      $$
      \omega = -KY\exp\left(-E_a/R_u T\right)
      $$
    
    + **In the scenario of supernova flame**, the internal energy is described by the equation of state of the relativistic degenerate electron gas[^2]:
      $$
      e{(\rho,T)}=\frac{3}{4{\rho}}(3\pi^2)^{\frac{1}{3}} \hbar c_1({\rho}N)^{\frac{4}{3}} + \frac{1}{2} N \frac{(3\pi^2)^{\frac{2}{3}}}{3 \hbar c_1}\left(\frac{1}{{\rho} N}\right)^{\frac{1}{3}} (k{T})^2
      $$
      With pressure being $p=\rho e/3$, one can get the form as: (constant $\gamma=4/3$)
      $$
      e=\frac{p}{\rho (4/3-1)} = \frac{p}{\rho (\gamma-1)}
      $$
      The reaction rate is adopted from the rate of Carbon-Carbon nuclear reaction[^3]:
      $$
      \omega = A Y^2 \exp({{- E_a}/{T_9^{\frac{1}{3}}}})
      $$
      where $T_9=T/10^9$ and $E_a=84.165$ for C12-C12 reaction.
    
      > $N$ is Carbon's electron density
      > $c_1$ is the light speed in vacuum
      > $\hbar$ is the reduced Plank constant
    
    

#### 2.2 Spherical flame

+ In spherical coordinates, the governing equations change to be:
  $$
  \begin{align*}
  \frac{\partial \boldsymbol U}{\partial t} + \frac{1}{r^2}\frac{\partial r^2 \boldsymbol  F(\boldsymbol U)}{\partial r} = \frac{1}{r^2}\frac{\partial r^2 \boldsymbol D(\boldsymbol U)}{\partial r} + \boldsymbol S
  \end{align*}
  $$
  or in the expanded form:
  $$
  \begin{align*}
  \frac{\partial \boldsymbol U}{\partial t} + \frac{\partial \boldsymbol  F(\boldsymbol U)}{\partial r} + \frac{2}{r}\boldsymbol  F(\boldsymbol U) = \frac{\boldsymbol D(\boldsymbol U)}{\partial r} + \frac{2}{r}\boldsymbol D(\boldsymbol U) + \boldsymbol S
  \end{align*}
  $$
  

  with the vectors being:
  $$
  \begin{align*}
  \boldsymbol U = \left[\begin{matrix} \rho \\ \rho u \\ \rho E \\ \rho Y\end{matrix}\right],
  \,\,\,\,
  \boldsymbol F(\boldsymbol U) = \left[\begin{matrix} \rho u \\ \rho u^2+p \\ \rho uH \\ \rho uY\end{matrix}\right],
  \,\,\,\,
  \boldsymbol D(\boldsymbol U) = \left[\begin{matrix} 0 \\ 0 \\  \lambda \frac{\partial T}{\partial r} + \rho QD \frac{\partial Y}{\partial r} \\  \rho D \frac{\partial Y}{\partial r} \end{matrix}\right],
  \,\,\,\,
  \boldsymbol S =\left[\begin{matrix} 0 \\ \frac{2p}{r} \\ 0 \\  \rho \omega \end{matrix}\right]
  \end{align*}
  $$
  

#### 2.3 Parameters for $H_2$-$O_2$ flame

+ Reference state ($H_2:O_2=2:1$)
  + pressure $p_0=1.013\times 10^5$ [$Pa$]
  + specific gas constant $R_u=692.25$ [$J/kg/K$]
  + temperature $T_0=300$ [$K$]
  + density $\rho_0=p_0/R_uT_0=0.4879$ [$kg/m^3$]
  + heat capacity ratio $\gamma=1.4$
  + velocity $u_0=\sqrt{R_uT_0}=455.7m/s$ ($S_L\approx10m/s$)
  + flame thickness $\delta_f\approx 5mm= 5\times 10^{-4} m$
  + time scale $t_0=\delta_f/u_0=1.1\times10^{-6}$ [$s$]
  + heat value $Q\approx 36 R_uT_0$ [$J/kg$] (300K->3076K)
  + thermal conductivity $\lambda\approx0.13$ [$W/m/K$]
  + thermal diffusivity $\alpha\approx 1\times10^{-4}$ [$m^2/s$]
  + Lewis number $Le\approx 0.3$
  + activation energy $Ea \approx 27 R_u T_0$
  + pre-exponential factor $A\approx 1\times10^7$ ?



### 3. Dimensionless Equations

#### 3.1 References variables

| Dimensionless variables | Reference values |
|--------------------|---------|
| $ \hat \rho = \frac{\rho}{\rho_0} $            | $ \rho_0 = 3.5\times 10^{10} kg/m^3 $ |
| $ \hat u  = \frac{u}{u_0} $                          | $ {u_0} = 7.871\times10^6 m/s,  S_L = 466m/s $ (laminar flame speed) |
| $ r = \frac{\hat r}{\hat {r_0}} $                                    | $ {r_0} = \delta_f = \frac{\lambda}{\rho_0 C_p S_t}=9\times10^{-4} m $ |
| $ \hat t = \frac{t}{t_0} $                               | $ t_0 = \frac{r_0}{u_0} = 1.143\times10^{-10} s$ |
| $ \hat T = \frac{T}{T_{b}-{T_0}} $ | $ {T_0}=1\times10^8 K; {T_{b}}=3.2\times 10^9 K $ |
| $ \hat p =  \frac{ p}{{p_0}} $                          | $ {p_0} = \rho_0u_0^2$ |
| $ \hat E = \frac{E}{p_0/\rho_0} = \frac{E}{u_0^2} $ | $e(\rho_0, T_b)= 6.1952\times 10^{13} J/kg $ |
| $ \hat Y = \frac{Y}{{Y_0}} $                            | $ {Y_0} = 1 $ |
| $ \hat\omega = \frac{ \omega}{ {\omega_0}} $ | $\rho_0{\omega_0}\delta_f=\rho_0Y_0 S_t$,  $ {\omega_0} = \frac{{Y_0}u_0}{r_0}\frac{S_t}{u_0} $ |
| $ \hat \alpha = \frac{\alpha}{u_0 r_0} $       | $\hat D= \frac{\hat \alpha}{Le}=\frac{2}{3}\frac{1}{Le}$ |
| $ \hat Q = \frac{Q}{u_0^2} $                  | $\hat Q=0.904$ |
| $ \hat \lambda = \lambda/(\frac{\rho_0 r_0u^3_0}{T_b-T_0})$ | $\hat \lambda=\frac{S_t}{u_0} \frac{C_p(T_b-T_0)}{u_0^2} = 0.14541$ |

> Laminar flame speed: $S_L=466 m/s$
>
> Turbulent flame speed: $S_t$, assumed to be equal to the speed of sound at burnt state

Therefore, the dimensionless equation can be obtained from: (use planar equation for demonstration)
$$
\begin{align*}
\frac{\partial \rho}{\partial t} + \frac{\partial }{\partial x}(\rho u) &= 0 \\
\frac{\partial (\rho u)}{\partial t} + \frac{\partial }{\partial x}(\rho u^2 +p) &= 0 \\
\frac{\partial \rho E}{\partial t} + \frac{\partial}{\partial x}(\rho uE+up) &= \lambda\frac{\partial^2 T}{\partial x^2} + \rho QD\frac{\partial^2 Y}{\partial x^2} \\
\frac{\partial \rho Y}{\partial t} + \frac{\partial }{\partial x}(\rho uY) &= \rho D\frac{\partial^2 Y}{\partial x^2} + \rho \omega
\end{align*}
$$

$\Rightarrow$
$$
\begin{align*}
\frac{\rho_0}{t_0}\left[\frac{\partial \hat \rho }{\partial \hat t} + \frac{\partial }{\partial \hat x}(\hat \rho \hat u)\right] &= 0 \\
\frac{\rho_0 u_0}{t_0} \left[\frac{\partial (\hat \rho \hat u)}{\partial \hat t} + \frac{\partial }{\partial \hat x}(\hat \rho \hat u^2 +\hat p)\right] &= 0 \\
\frac{\rho_0 u_0^2}{t_0}\left[\frac{\partial \hat \rho \hat E}{\partial \hat t} + \frac{\partial}{\partial \hat x}(\hat \rho \hat u\hat E+\hat u\hat p)\right] &= 
\frac{\rho_0u_0^3r_0(T_b-T_0)}{r_0^2(T_b-T_0)} \left[\hat \lambda\frac{\partial^2 \hat T}{\partial \hat x^2}\right] + \frac{\rho_0u_0^2u_0r_0Y_0}{r_0^2}\left[\hat \rho \hat Q \hat D\frac{\partial^2 \hat Y}{\partial \hat x^2}\right] \\
\frac{\rho_0Y_0}{t_0} \left[\frac{\partial \hat \rho \hat Y}{\partial \hat t} + \frac{\partial }{\partial \hat x}(\hat \rho \hat u \hat Y)\right] &= \frac{\rho_0 u_0r_0Y_0}{r_0^2}\left[\hat \rho \hat D\frac{\partial^2 \hat Y}{\partial \hat x^2}\right] + \frac{\rho_0 u_0 Y_0}{r_0}\frac{S_t}{u_0}[\hat \rho \hat \omega]
\end{align*}
$$

$\Rightarrow$
$$
\begin{align*}
\frac{\partial \hat \rho }{\partial \hat t} + \frac{\partial }{\partial \hat x}(\hat \rho \hat u) &= 0 \\
\frac{\partial (\hat \rho \hat u)}{\partial \hat t} + \frac{\partial }{\partial \hat x}(\hat \rho \hat u^2 +\hat p) &= 0 \\
\frac{\partial \hat \rho \hat E}{\partial \hat t} + \frac{\partial}{\partial \hat x}(\hat \rho \hat u \hat H) &= 
\hat \lambda\frac{\partial^2 \hat T}{\partial \hat x^2} + \hat \rho \hat Q \hat D\frac{\partial^2 \hat Y}{\partial \hat x^2} \\
\frac{\partial \hat \rho \hat Y}{\partial \hat t} + \frac{\partial }{\partial \hat x}(\hat \rho \hat u \hat Y) &= \hat \rho \hat D\frac{\partial^2 \hat Y}{\partial \hat x^2} + \frac{S_t}{u_0}\hat \rho \hat \omega
\end{align*}
$$

#### 3.2 Dimensionless parameters

+ Reference velocity $u_0$:

  For simplification, the largest internal energy is normalized to be $\hat e(\rho,T_b) = e({\rho_0},{T_b})/ u_0^2 = 1$. Therefore, one has
  $$
  u_0 = \sqrt{e(\rho_0,T_b)} = \sqrt{6.1952\times10^{13}} = 7.871 \times 10^6 m/s
  $$

+ Reference length scale $r_0=9\times10^{-4}m$

+ Reference time scale $t_0=r_0/u_0=1.143\times10^{-10}s$

+ Speed of sound at the burnt state $c= \sqrt{\gamma p/\rho} =\sqrt{\gamma(\gamma-1)e(\rho_0,T_b)} = 5.247\times10^6 m/s$

+ Heat value $\hat Q=Q/u_0^2$,  with data from Fowler 1975, for $C^{12}+C^{12} \to Mg^{24}$ , $Q=13.931MeV$, thus
  $$
  \begin{align*}
  \hat Q 
  &= \frac{13.931MeV \times (1.6022\times10^{-13} J/MeV) \times (6.022\times 10^{23}  mol^{-1})/(0.024kg/mol)}{u_0^2}\\
  &= \frac{5.6\times10^{13}J/kg}{6.1952\times10^{13}J/kg} = 0.904
  \end{align*}
  $$

+ The speed ratio $S_t/u_0 = 2/3$
  
+ Thermal conductivity $\hat \lambda$ is related to $C_p$, while the heat release is related to internal energy:
  $$
  \tilde C_p(T_b-T_0) = \int_{T_0}^{T_b} C_p dT = e(\rho_0,T_b) - e(\rho_0, T_0)
  $$

  $$
  \hat \lambda = \frac{S_t}{u_0} \frac{\tilde C_p(T_b-T_0)}{u_0^2} = 0.14541
  $$

+ The diffusivity $D$ is described by Lewis number $D = \alpha / Le = \lambda/\rho C_pLe$, therefore
  $$
  \hat D = \frac{\lambda}{\rho C_pLe} \frac{1}{u_0 r_0} = \frac{1}{Le}\frac{S_t }{u_0} \frac{\delta_f}{r_0} = \frac{2}{3}\frac{1}{Le}
  $$

+ The total energy form keeps  the same in dimensionless form:
  $$
  \hat E = \hat e + \frac{1}{2}\hat u^2 + \hat Q\hat Y, \, \hat e = \frac{\hat p}{\hat \rho (\gamma-1)}
  $$
  with the equation of state as:
  $$
  \begin{align*}
  p &= \rho e(\gamma-1) 
  = \frac{3}{4}(\gamma-1)(3\pi^2)^{\frac{1}{3}} \hbar c_1({\rho}N)^{\frac{4}{3}} + \frac{1}{2} (\gamma-1)N \frac{(3\pi^2)^{\frac{2}{3}}}{3 \hbar c_1}\left(\frac{{\rho}^2}{ N}\right)^{\frac{1}{3}} (k{T})^2 \\
  
  \hat p
  & = \frac{3}{4}\frac{1}{\rho_0 u_0^2}(\gamma-1)(3\pi^2)^{\frac{1}{3}} \hbar c_1({\rho_0}N)^{\frac{4}{3}} \hat \rho^{\frac{4}{3}} + \frac{1}{2} (\gamma-1)N \frac{(3\pi^2)^{\frac{2}{3}}}{3 \hbar c_1}\left(\frac{{\rho_0}^2}{ N}\right)^{\frac{1}{3}} k^2(T_b-T_0)^2 \hat \rho^{\frac{2}{3}} \hat T^2 \\
  &= C_1\hat\rho^{4/3}+C_2\hat\rho^{2/3}\hat T^2
  \end{align*}
  $$
  with $C_1=2.60559\times10^{-1}$, $C_2=6.82973\times10^{-2}$

+ The source term equation:
  $$
  \omega_0 = \frac{{Y_0}u_0}{r_0}\frac{S_t}{u_0} \approx A Y_0^2 \exp({{- E_a}/{T_{b9}^{\frac{1}{3}}}})  \\
  \hat \omega = \frac{\omega}{\omega_0}= \hat Y^2 \exp(-Ea/T_9^{\frac{1}{3}}+Ea/T_{9b}^{\frac{1}{3}}) = \hat Y^2\exp\left( -Ea\left(\frac{T_b-T_0}{10^9}\right)^{-\frac{1}{3}}\left(\hat T^{-\frac{1}{3}}-\hat T_b^{-\frac{1}{3}}\right)\right)
  $$
  to reduce the term $S_t/u_0$ in previous equation along with $\hat \omega$, denoting
  $$
  \begin{align*}
  \widehat{Ea} &= -Ea\left(\frac{T_b-T_0}{10^9}\right)^{-\frac{1}{3}} = 57.7224 \\
  \hat A &= \frac{S_t}{u_0}\exp\left(Ea\left(\frac{T_b-T_0}{10^9}\right)^{-\frac{1}{3}}\hat T_b^{-\frac{1}{3}}\right) = 4.25133 \times10^{24}
  \end{align*}
  $$
  then the source term keeps the same in dimensionless form:
  $$
  \hat \omega = \hat A \hat Y^2 \exp(-\widehat{Ea}/\hat T^{\frac{1}{3}} )
  $$

+ And the final dimensionless equations are:
  $$
  \begin{align*}
  &\frac{\partial \hat \rho }{\partial \hat t} + \frac{\partial }{\partial \hat x}(\hat \rho \hat u) = 0 \\
  &\frac{\partial (\hat \rho \hat u)}{\partial \hat t} + \frac{\partial }{\partial \hat x}(\hat \rho \hat u^2 +\hat p) = 0 \\
  &\frac{\partial \hat \rho \hat E}{\partial \hat t} + \frac{\partial}{\partial \hat x}(\hat \rho \hat u \hat H) = \hat \lambda\frac{\partial^2 \hat T}{\partial \hat x^2} + \hat \rho \hat Q \hat D\frac{\partial^2 \hat Y}{\partial \hat x^2} \\
  & \frac{\partial \hat \rho \hat Y}{\partial \hat t} + \frac{\partial }{\partial \hat x}(\hat \rho \hat u \hat Y) = \hat \rho \hat D\frac{\partial^2 \hat Y}{\partial \hat x^2} + \hat \rho \hat \omega \\
  &\hat E = \hat e + \frac{1}{2}\hat u^2 + \hat Q\hat Y \\
  &\hat e = \frac{\hat p}{\hat \rho (\gamma-1)} \\
  &\hat p = C_1\hat\rho^{4/3}+C_2\hat\rho^{2/3}\hat T^2 \\
  &\hat \omega = \hat A \hat Y^2 \exp(-\widehat{Ea}/\hat T^{\frac{1}{3}} )
  \end{align*}
  $$



### 4. Numerical methods

+ Discrete methods:

  + Time evolution: 3 order TVD `Runge-Kutta`
    $$
    \begin{align*}
    U^{(1)} &= U_n + L(U_n) \Delta t \\
    U^{(2)} &= \frac{3}{4}U_n + \frac{1}{4}(U^{(1)} + L(U^{(1)}) \Delta t) \\
    U_{n+1} &= \frac{1}{3}U_n + \frac{2}{3}(U^{(2)} + L( U^{(2)})\Delta t) \end{align*}
    $$

  + Convection term:

    `Roe` method: solve convection flux by eigen vector

  + Diffusion term: 7 order central difference
    $$
    \frac{\partial m_i}{\partial x} = \frac{1}{\Delta x}\left(\frac{1}{60}m_{i+3}-\frac{9}{60}m_{i+2}+\frac{45}{60}m_{i+1}-\frac{45}{60}m_{i-1}+\frac{9}{60}m_{i-2}-\frac{1}{60}m_{i-3} \right)
    $$

+ Initial conditions:

  + Pressure: whole field the same pressure $\hat P(\rho_0, T_0)/\hat P_0$
  + Velocity: $u=0$
  + Concentration: $Y_0=1$
  + Temperature: $T_0 \to T_b$, $\tanh(x)$ profile for about 1 flame thickness area

+ Boundary conditions: (cartesian coordinates)

  + Inner boundary conditions: `opening / freeflow / inletOutlet`
  + Outlet boundary conditions: `opening / freeflow / inletOutlet`

+ Boundary conditions: (spherical coordinates)

  + 
  



### 5. Results





### Appendix

#### A1. Speed of sound

+ Pressure from equation of state:

  $$
  \begin{align*}
  p = \rho e(\gamma-1) 
  &= \frac{3}{4}(\gamma-1)(3\pi^2)^{\frac{1}{3}} \hbar c_1({\rho}N)^{\frac{4}{3}} + \frac{1}{2} (\gamma-1)N \frac{(3\pi^2)^{\frac{2}{3}}}{3 \hbar c_1}\left(\frac{{\rho}^2}{ N}\right)^{\frac{1}{3}} (k{T})^2 \\
  &= C_1 \rho^{\frac{4}{3}} + C_2 \rho^{\frac{2}{3}}
  \end{align*}
  $$

+ The theoretical speed of sound is obtained by:
  $$
  c 
  = \sqrt{ \left(\frac{\partial p}{\partial \rho} \right)_s}
  = \sqrt{ \frac{4}{3}C_1\rho^{\frac{1}{3}} + \frac{2}{3}C_2\rho^{-\frac{1}{3}} } \approx \sqrt{\gamma \frac{p}{\rho}}
  $$

####  A2. Eigenvector of Jacobian


+ Eigen vector for `Roe` method:
  $ U = [\rho, \rho u, \rho E, \rho Y]  = [m_1, m_2, m_3, m_4] $
  $ F = [\rho u, \rho u^2 + p, \rho u E + u p, \rho u Y] = [m_2, m_2^2/m_1+p，(m_3 + p) m_2 / m_1, m_4 m_2/m_1 ] $

  where $p=(\gamma-1)\rho(E-\frac{1}{2}u^2-QY)$

  + Solve by `Matlab` symbolic calculation:
  
    ```matlab
    % definition
    syms m1 m2 m3 m4 real
    syms Q gamma real
    syms m1 m3 m4 gamma Q positive
    
    % declaration
    rho = m1;
    u = m2/m1;
    E = m3/m1;
    Y = m4/m1;
    e = E-1/2*u*u-Q*Y
    p = (gamma-1)*rho*e;
    Um = [m1 m2 m3 m4];
    Fm = [rho*u, rho*u*u+p; rho*u*E+u*p; rho*u*Y];
    
    % solve
    J = jacobian(Fm,Um);
    [RM,RD] = eig(J);
    [LM,LD] = eig(J');
    ```
  
  + The Jacobian matrix is:
    $$
    \boldsymbol J = \left[\begin{matrix} 
    0          & 1        & 0        & 0    \\
    -u^2 + K_1  & 2u +K_2  & K_3      & K_4  \\
    -uH + uK_1 & H+uK_2   & u(1+K_3) & uK_4 \\
    -uY        & Y        & 0        & u
    \end{matrix}\right]
    $$
    where $K_1 = (\gamma-1)u^2/2,\ K_2 = (\gamma-1)(-u),\ K_3=(\gamma-1),\ K_4=(\gamma-1)(-Q) $
    
  + Eigenvalue $\Lambda$ is:
    $$
    \begin{aligned}
    \boldsymbol \Lambda =
    \left[\begin{matrix}
    u - c       \\
    & u         \\
    & & u       \\
    & & & u + c
    \end{matrix}\right]
    \end{aligned}
    $$
    where $c$ is the numerical speed of sound:
    $$
    c = \frac{2}{3}\sqrt{E - \frac{1}{2}u^2 - QY} = \frac{2}{3}\sqrt{\frac{p}{\rho}\frac{1}{\gamma-1}} = \sqrt{\gamma p/\rho}
    $$
    
  + Right eigen matrix $\boldsymbol R$ is: (each column is an eigenvector)
    $$
    \boldsymbol R = 
    \left[\begin{matrix}
    1    & 1    & 0 & 1       \\
    u-c  & u    & 0 & u+c     \\
    H-uc & \color{red}{u^2/2} & Q  & H + uc  \\
    Y    & 0    & 1 & Y
    \end{matrix}\right]
    $$
  
  + Left eigen matrix $\boldsymbol L$ is: (each row is an eigenvector)
    $$
    \boldsymbol L = 
    \left[\begin{matrix}
    -\frac{u}{2c}-\frac{K_1}{2c^2} & \frac{1}{2c}+\frac{uK_3}{2c^2} & -\frac{K_3}{2c^2} & -\frac{K_4}{2c^2} \\
    -\frac{1}{2} + \frac{K_1}{2c^2}  & -\frac{uK_3}{2c^2} & \frac{K_3}{2c^2} & \frac{K_4}{2c^2} \\
    -Y & 0 & 0 & 1  \\
    -\frac{u}{2c}+\frac{K_1}{2c^2} & \frac{1}{2c}-\frac{uK_3}{2c^2} & \frac{K_3}{2c^2} & \frac{K_4}{2c^2}
    \end{matrix}\right]
    $$
  
  + To meet the eigen-decomposition criteria and to be used in the eigen-space projection, re-projection step, one can get the eigen-space as: (let $\gamma_1\equiv\gamma-1$, then $K_1=\gamma_1u^2/2$, $K_3=\gamma_1$, $K_4=-\gamma_1Q$)
    $$
    \boldsymbol{R_F} = 
    \left[\begin{matrix}
    \frac{1}{c}    & \frac{1}{c}    & 0           & \frac{1}{c}      \\
    \frac{u}{c}-1  & \frac{u}{c}    & 0           & \frac{u}{c}+1    \\
    \frac{H}{c}-u  & \frac{u^2}{2c} & \frac{Q}{c} & \frac{H}{c} + u  \\
    \frac{Y}{c}    & 0              & \frac{1}{c} & \frac{Y}{c}
    \end{matrix}\right]
    ,\
    \boldsymbol{R_F^{-1}} = 
    \left[\begin{matrix}
    \frac{u}{2}+\frac{\gamma_1u^2}{4c} & -\frac{1}{2}-\frac{\gamma_1u}{2c} & \frac{\gamma_1}{2c} & -\frac{\gamma_1Q}{2c} \\
    c - \frac{\gamma_1u^2}{2c}         & \frac{\gamma_1u}{c}               & -\frac{\gamma_1}{c} &  \frac{\gamma_1Q}{c}  \\
    -\frac{\gamma_1 u^2Y}{2c}          & \frac{\gamma_1 uY}{c}             & -\frac{\gamma_1Y}{c}&c+\frac{\gamma_1QY}{c} \\
    -\frac{u}{2}+\frac{\gamma_1u^2}{4c}& \frac{1}{2}-\frac{\gamma_1 u}{2c} & \frac{\gamma_1}{2c} & -\frac{\gamma_1Q}{2c}
    \end{matrix}\right]^T
    $$
    which satisfies $\boldsymbol{R_F} \boldsymbol{R_F}^{-1} = \boldsymbol{I}$ and $\boldsymbol{R_F}\boldsymbol{\Lambda}\boldsymbol{R_F}^{-1} = J$.
  
  + This has been validated in `Matlab` code:
  
    ```matlab
    syms u real
    syms rho E Y real positive
    syms Q gamma real positive
    
    e = (E - 1/2*u*u - Q*Y);
    p = (gamma-1)*rho*e;
    c = sqrt(gamma*p/rho);
    H = E + p/rho;
    
    K1 = (gamma-1)*u*u/2;
    K2 = (gamma-1)*(-u);
    K3 = (gamma-1);
    K4 = (gamma-1)*(-Q);
    
    J = [
      [        0,      1,         0,    0];
      [  -u*u+K1, 2*u+K2,        K3,   K4]
      [-u*H+u*K1, H+u*K2,  u*(1+K3), u*K4]
      [     -u*Y,      Y,         0,    u];
    ];
    
    D = [
      [u-c, 0, 0,   0];
      [  0, u, 0,   0];
      [  0, 0, u,   0];
      [  0, 0, 0, u+c];
    ];
    
    RF = [
      [1/c,       1/c,      0,     1/c];
      [u/c-1,     u/c,      0,   u/c+1];
      [H/c-u, u*u/2/c,    Q/c,   H/c+u];
      [Y/c,         0,    1/c,     Y/c];
    ];
    
    g1 = gamma-1;
    c2 = 2*c;
    c4 = 4*c;
    LF = [
      [ u/2+g1*u*u/c4,   -1/2-g1*u/c2,   g1/c2,    -g1*Q/c2];
      [   c-g1*u*u/c2,        g1*u/c,   -g1/c,      g1*Q/c ];
      [  -Y*g1*u*u/c2,      Y*g1*u/c, -Y*g1/c,  c+Y*g1*Q/c ];
      [-u/2+g1*u*u/c4,    1/2-g1*u/c2,   g1/c2,    -g1*Q/c2];
    ]';
    
    simplify(RF*LF')        % output is [ 1, 0, 0, 0]
                            %           [ 0, 1, 0, 0]
                            %           [ 0, 0, 1, 0]
                            %           [ 0, 0, 0, 1]
                            
    simplify(RF*D*LF' - J)  % output is [ 0, 0, 0, 0]
                            %           [ 0, 0, 0, 0]
                            %           [ 0, 0, 0, 0]
                            %           [ 0, 0, 0, 0]
    ```
  



### Reference

[^1]: Wheeler J C, Harkness R P, Rep. Prog. Phys. 1990, 53:1467-1557
[^2]: G. Xing, Y. Zhao, et al., Astro. J., 2017, 841:21
[^3]: Fowler, W. A., Caughlan, G. R., & Zimmerman, B. A. 1975, ARA&A, 13, 69

[3] https://en.wikipedia.org/wiki/Divergence#Spherical_coordinates<br>[4] Landau L D , Lifshitz E M . Statistical Physics, Part 1[J]. Physics Today, 1980.<br>[5] Woosely. 2011. FLAMES IN TYPE Ia SUPERNOVA: DEFLAGRATION–DETONATION TRANSITION IN THE OXYGEN-BURNING FLAME<br>[6] http://www.astrophysicsspectator.org/topics/stars/FusionCarbonOxygen.html