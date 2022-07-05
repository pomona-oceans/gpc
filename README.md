# gpc

The purpose of this project is to explore the role of global phosphorous cycling in the unusual geochemical landscapes of the late Neoproterozoic. We aim to account for a variety of factors including biogeochemistry, erosion, and sedimentary deposition within their tectonic context.

Our strategy is to construct a simple box model representing components and interrelationships of the global phosphorous reservoir. This model follows, e.g., chapter 5 of Chameides & Perdue's _Biogeochemical Cycles._ The basic construction is a general linear system of the form:

$$\mathbf{y}'(t) = A(t)\mathbf{y}(t) + \mathbf{b}(t).$$

For a system with $n$ boxes, $\mathbf{y}(t)$ and $\mathbf{y}'(t)$ are vectors with $n$ entries. $A(t)$ is an $n\times n$ matrix describing the linear interdependencies of each box, and $\mathbf{b}(t)$ is an $n-$vector meant to capture external "forcing."

The simplest starting case is unforced, which means $\mathbf{b}(t)$ = $\overrightarrow{0}$. Additionally, we assume the linear dependencies in $A(t)$ are constant, i.e., $A(t)=A$. Thus, the system is given by:

$$\mathbf{y}'(t) = A\mathbf{y}(t).$$

Given initial conditions $\mathbf{y}(0) = \mathbf{y}_0$, the solution to this system is 

$$\mathbf{y}(t) = e^{At}\cdot\mathbf{y}_0.$$
## Derivation

A three-box model captures most of the complexity of a general $n$-box model, so we'll use that to explore the procedure. Following Chameides & Perdue (1997), I use the capital letter $F$ to indicate fluxes of phosphate between reservoirs in teragrams per year.

Note on ordering conventions: This is very important because each flux represents an _ordered_ relationship between two of the three (or $n$) boxes. That is, the flux from A to B is not the same as from B to A. Chameides & Perdue (1997) use $F_{i\to j}$ as the flux from box $i$ to box $j$. If each row corresponds to a scalar differential equation, this means $F_{i\to j}$ goes in row $j$ and column $i$. However, the convention in linear algebra (and python) is to express entries of a matrix as an ordered pair $F_{ij}$ where $i$ is the row and $j$ is the column. This is $F[i,j]$ in python. To avoid confusion, I have included some reminders to always define entries in order [TO, FROM] rather than think about rows and columns each time. To adapt Chameides & Perdue, think $F_{i\leftarrow j}$.

In a general system of three boxes, any ordered pair of two different boxes may be assigned a unique entry. The important point is that diagonal entries (i.e., $i=j$) are not interpretable in our model. As a matrix, this looks like:

$$F = \begin{bmatrix}
  & f_{0\leftarrow1} &   f_{0\leftarrow2}\\
    f_{1\leftarrow0} & & f_{1\leftarrow2}\\
    f_{2\leftarrow0} &   f_{2\leftarrow1}
\end{bmatrix}
(\text{Tg P yr}^{-1})$$

Notice that numbering begins at 0, not 1, which is the convention in most programming languages. So far, we haven't really done anything besides arrange the fluxes of the model in a convenient way. And the numbers we put in here are most likely derived from some combination of prior modeling and observation to reflect a particular moment in time, such as modern day. These actual fluxes are certainly dependent on the actual size of modern day reservoirs, for example.

We want to take these specific numbers and transform them into something more general which we can modify and experiment with. One of the simplest ways to generalize this system is to assume each flux is directly proportional to the reservoir from which it comes.

To do that, notice that each entry represents the flux from the reservoir corresponding to the column it's in. That means we could represent each flux entry (with units Tg P /yr) as a product of a reservoir (Tg P) and a linear _rate constant_ (/yr). It's convenient that we have fairly good estimates for the size of global phospate reservoirs. For the purposes of this illustration, let's arrange these reservoirs in a simple array:

$$\mathbf{p} = \begin{bmatrix}
  p_0 \\ p_1 \\ p_2
\end{bmatrix}(\text{Tg P})$$

That means our rate constant matrix will be constructed as follows:
 $$K = \begin{bmatrix}
& f_{0\leftarrow1} / p_1 &   f_{0\leftarrow2} / p_2\\
  f_{1\leftarrow0} / p_0 & & f_{1\leftarrow2} / p_2\\
  f_{2\leftarrow0} / p_0 &   f_{2\leftarrow1} / p_1
\end{bmatrix} = \begin{bmatrix}
  & k_{0\leftarrow1} &   k_{0\leftarrow2}\\
    k_{1\leftarrow0} & & k_{1\leftarrow2}\\
    k_{2\leftarrow0} &   k_{2\leftarrow1}
\end{bmatrix} (\text{yr}^{-1}).$$

The point of doing this factorization is that we might want to experimentally manipulate $K$ and $\mathbf{p}$ as a function of time $t$, and see how they affect $F$ and subsequenctly the model behavior. Remember, $K$ and $\mathbf{p}$ have physical meaning: one being the magnitude of the phosphate reservoirs, and the other being the _sensetivity_ of each flux to reservoir magnitudes. We certainly expect the reservoir magnitudes to change over time as phosphate fluxes throughout the system. But we might also expect that global climate or tectonic conditions might change the "sensitivity" as well. For example, even holding all reservoir magnitudes constant, a hotter climate will tend to promote more rapid weathering and erosion, so $K$ could be a function of $t$.

There's just one glaring problem here. We first defined our matrices in such a way that the diagonal terms are zero, since we can't really interpret a flux from a reservoir to itself. But remember, we can't really make any time-dependent predictions using just $F$ or even with $K\mathbf{p}$. We want an expression for how each reservoir (magnitude) evolves over time, which is to say $d\mathbf{p}/dt$ or $\mathbf{p}'$. We would love something like

$$\mathbf{p}'=K\mathbf{p}$$

which we (the computer) can solve numerically or analytically. To see the problem with our construction so far, we can write out this system explicitly for each reservoir:
$$
  \mathbf{p}'_0=
  \mathbf{p}_1k_{0\leftarrow1} + \mathbf{p}_2k_{0\leftarrow2}\\
  \mathbf{p}'_1=
  \mathbf{p}_0k_{1\leftarrow0}  + \mathbf{p}_2k_{1\leftarrow2}\\
  \mathbf{p}'_2=  
  \mathbf{p}_0k_{2\leftarrow0} + \mathbf{p}_1k_{2\leftarrow1}
$$
You should check that this system is exactly equivalent to $\mathbf{p}'=K\mathbf{p}$ as $K$ has been defined so far. But there's something missing. We need to account for how each reservoir affects its own rate of change. Since each flux is proportional to the size of its source, we should expect that as a reservoir gets larger, its rate of change becomes more negative, since the rate of outflux is proportional to magnitude. For each flux, we have already accounted for the positie influence it has on its sink. Now we just need to match these to the negative influence it has on its source:
$$\begin{gather*}
  \mathbf{p}'_0=
  \mathbf{p}_0(-k_{1\leftarrow0}-k_{2\leftarrow0}) + \mathbf{p}_1k_{0\leftarrow1} + \mathbf{p}_2k_{0\leftarrow2}\\
  \mathbf{p}'_1=
  \mathbf{p}_0k_{1\leftarrow0} + \mathbf{p}_1(-k_{0\leftarrow1}-k_{2\leftarrow1}) + \mathbf{p}_2k_{1\leftarrow2}\\
  \mathbf{p}'_2=  
  \mathbf{p}_0k_{2\leftarrow0} + \mathbf{p}_1k_{2\leftarrow1} + \mathbf{p}_2(-k_{0\leftarrow2}-k_{1\leftarrow2})
\end{gather*}$$

That means we need to adjust the matrix $K$ (and thus $F$) in a very particular way. To see how, let's express this correct system in terms of a matrix vector product:

$$\begin{bmatrix}
  p_0 \\ p_1 \\ p_2
\end{bmatrix}'=\begin{bmatrix}
  -k_{1\leftarrow0}-k_{2\leftarrow0} & k_{0\leftarrow1} &   k_{0\leftarrow2}\\
  k_{1\leftarrow0} & -k_{0\leftarrow1}-k_{2\leftarrow1} & k_{1\leftarrow2}\\
  k_{2\leftarrow0} & k_{2\leftarrow1} & -k_{0\leftarrow2}-k_{1\leftarrow2}
\end{bmatrix}\begin{bmatrix}
  p_0 \\ p_1 \\ p_2
\end{bmatrix}.$$

You may notice a few convenient facts. First, each $k_{i\leftarrow j}$ term is written exactly twice, once negative for its source and once positive for its sink. Second, these entries are arranged in such a way that each column of $K$ sums to zero. That means there's a convenient and easy way to transform from our first version (the one with zeroes on the diagonal) and the correct version here. If we start with

$$\begin{bmatrix}
  & k_{0\leftarrow1} &   k_{0\leftarrow2}\\
    k_{1\leftarrow0} & & k_{1\leftarrow2}\\
    k_{2\leftarrow0} &   k_{2\leftarrow1}
\end{bmatrix},$$

we can imagine a heavy weight sitting on top, squishing the whole matrix into a flat row vector with entries corresponding to the sum of each column:

$$\begin{bmatrix}
(k_{1\leftarrow0} + k_{2\leftarrow0}) & (k_{0\leftarrow1} + k_{2\leftarrow1}) & (k_{0\leftarrow2} + k_{1\leftarrow2})
\end{bmatrix}.$$

Then we want to tip this vector up 45 degrees to the right to form a diagonal matrix like so:

$$\begin{bmatrix}
(k_{1\leftarrow0} + k_{2\leftarrow0}) & 0 & 0\\ 0 & (k_{0\leftarrow1} + k_{2\leftarrow1}) & 0\\ 0 & 0 & (k_{0\leftarrow2} + k_{1\leftarrow2})
\end{bmatrix}$$

with zeroes in all the other entries. Now we have just the right missing puzzle piece to fill in the gaps of our first attempt. We just have to make sure to subtract, rather than add, to match the signs we want.

And that's it! The net result is to take a matrix and adjust its diagonal entries to ensure each column sums to zero. In this case, we started with a very specific matrix—one whose diagonal entries are zero. I call this a "hollow" matrix and have defined the whole procedure as "fill_hollow." But the procedure works in general. As a check on your understanding, it would be good to predict and confirm what the effect is of completing this operation twice on the same matrix; that is, to fill a hollow matrix and then fill the result again.

In python, the squishing is called ".sum(axis=0)" and the tipping is called "np.diag." Alternatively, ".sum(axis=1)" would do the same thing for rows instead of columns, squishing from the left.
