# Weighted Structured Lasso SVM (WSLSVM)
SVM + Lasso + Structured Lasso 

hinge loss function: $$L(\beta_0,\beta; X, Y) = \frac{1}{n}\sum_{i=1}^n[1-y_{i}(\beta_0 + x_i^\top\beta )]_{+},$$

Lasso: $$\frac{\lambda\alpha}{n}\Vert D_{\omega_1}\beta \Vert_1,$$

Structured Lasso : $$\frac{\lambda(1-\alpha)}{n}\Vert D_{\omega_2}\mathbf{G}\beta \Vert_1,$$
where $X=(x_1,\cdots,x_n)^\top$, $Y=(y_1,\cdots,y_n)^\top$, $\beta = (\beta_1,\cdots,\beta_p)^{\top}$, $\mathbf{1}$ is an $n$-dimensional column of ones. The weight vectors $\omega_1$ and $\omega_2$ are associated with $\beta$ and $\mathbf{G}$, respectively, and $D_\omega$ denotes a diagonal matrix with elements from vector $\omega$.
