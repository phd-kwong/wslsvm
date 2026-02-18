# WSLSVM
$$(\hat{\beta}_0,\hat{\beta}) =
\argmin \limits_{(\beta_0,\beta)\in \mathbb{R}^{p+1}} L(\beta_0,\beta; X, Y) + \frac{\lambda\alpha}{n}\Vert D_{\omega_1}\beta \Vert_1+ \frac{\lambda(1-\alpha)}{n}\Vert D_{\omega_2}\mathbf{G}\beta \Vert_1,$$

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
