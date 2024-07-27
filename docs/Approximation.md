# Notes

An ellipsoid $\mathcal{E}_{Q,c}$ in $\mathbb{R}^d$ is specified by a $d\times d$ symmetric positive definite matrix $Q$ and a center $c \in \mathbb{R}$ and is defined as

$$
\mathcal{E}_{Q,c} = \left\{ x \in \mathbb{R}^d : (x-c)^T Q(x-c) \leq 1 \right\}
$$

The volume is given by
$$
\text{vol }\mathcal{E}_{Q,c} =  \eta \text{det} Q^{-1/2}
$$
where $\eta$ is the volum of the unit ball in $\mathbb{R}^d$.

It was proved (John) that




## For points

Given $\epsilon>0$, an ellipsoid $E_{Q,c}$ is said to be a $(1+\epsilon)$-approximation to $\rm MVEE (\mathcal{S})$ if
$$
E_{Q,c} \supseteq \mathcal{S}, \quad \dfrac{\text{vol } E_{Q,c}}{\text{vol MVEE}(\mathcal{S})} \leq (1+\epsilon) 
$$


## The KY algorithm for points

Two phases.
1. Initial volume approximation. Runs the Betke and Henk algorithm to collect a set $\mathcal{X}_0$ of $2d$ points.
2. These set of $2d$ points are refined by using Khachiyan's algoriothm to produce a $(1+\epsilon)$-approximation to $\rm MVEE (\mathcal{S})$.

### Betke and Henk algorithm

Gives a very rough approximation oto the volume of the $\mathcal{C H}(\mathcal{S})$ (the convex hull) but at the same time it produces the approximation by using at most a set of $2d$ vertices of the $\mathcal{CH}(\mathcal{S})$.

1. Choose a random direction $b_1$ in $R^d$ and find supposting 



The resulting set $\mathcal{X}_0$ satisfies that
$$
\text{vol conv}(\mathcal{S}) \leq d!\text{vol conv}(\mathcal{X}_0) 
$$