# Weighted Structured Lasso SVM (WSLSVM)
SVM + Lasso + Structured Lasso 

hinge loss function: $$L(\beta_0,\beta; X, Y) = \frac{1}{n}\sum_{i=1}^n[1-y_{i}(\beta_0 + x_i^\top\beta )]_{+},$$

Lasso: $$\frac{\lambda\alpha}{n}\Vert D_{\omega_1}\beta \Vert_1,$$

Structured Lasso : $$\frac{\lambda(1-\alpha)}{n}\Vert D_{\omega_2}\mathbf{G}\beta \Vert_1,$$
where $X=(x_1,\cdots,x_n)^\top$, $Y=(y_1,\cdots,y_n)^\top$, $\beta = (\beta_1,\cdots,\beta_p)^{\top}$, $\mathbf{1}$ is an $n$-dimensional column of ones. The weight vectors $\omega_1$ and $\omega_2$ are associated with $\beta$ and $\mathbf{G}$, respectively, and $D_\omega$ denotes a diagonal matrix with elements from vector $\omega$.

# ADMM
\begin{algorithm}
\caption{Alternating Direction Method of Multipliers (ADMM) for WSLSVM}
\label{alg:ADMM-WSLSVM}
\textbf{Input:} training data $(x_i,y_i)_{i=1}^{n}$, structural matrix $\mathbf{G}$, parameters $\lambda$, $\alpha$, $\nu$, $\rho$, tolerance $\epsilon^{abs}>0$, $\epsilon^{rel}>0$, maximum iteration number $M$
\begin{algorithmic}[1]
\State Initialize $z^k$, $b^k$, $\gamma^k$, $u^k$, $s^k$, $t^k$, $k=0$
\State Calculate $\omega_j=\tfrac{\Vert \widetilde{X}_{\cdot j}\Vert_2}{\sqrt{n}}$ for $j = 1, \ldots, p+m$, where $\omega = (\omega_1^\top, \omega_2^{\top})^\top$.
\State \textbf{Repeat}
\State $Z^{k+1}=\tau_{\frac{1}{n\nu}}( z^k-u^k)$.
\State $\beta^{k+1}=\mathcal{S}_{\tfrac{\lambda\alpha\omega_{1}}{n\nu}}(b^k-s^k)$s.
\State $\Gamma^{k+1}=\mathcal{S}_{\tfrac{\lambda(1-\alpha)\omega_{2}}{n\nu}}(\gamma^k-t^k)$.
\State $b^{k+1}=\Big((YX)^\top YX + \mathbf{G}^\top\mathbf{G}+I_p\Big)^{-1}E$.
\State $z^{k+1}=\mathbf{1}-YXb^{k+1}$.
\State $\gamma^{k+1}=\mathbf{G}b^{k+1}$.
\State $u^{k+1}=u^k+(Z^{k+1}-z^{k+1})$.
\State $s^{k+1}=s^k+(\beta^{k+1}-b^{k+1})$.
\State $t^{k+1}=t^k+(\Gamma^{k+1}-\gamma^{k+1})$.
\State $k\leftarrow k+1$.
\State \textbf{Until} stopping criterion is satisfied.
\end{algorithmic}
\textbf{Output:} the solution $\hat{\beta}=b^{k}$.
\end{algorithm}
