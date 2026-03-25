#import "@preview/subpar:0.2.2"
#import "lib.typ": project
#import "@preview/dashy-todo:0.1.3": todo
#import "@preview/wrap-it:0.1.1": *
#set text(font: "xits")
#set cite(style: "ieee")
#show: project.with(
  title: "Homework",
  subtitle: "Perturbation Methods MATH 2015-1",
  authors: (
    "Loïc Delbarre",
    "s215072"
  ),
  school-logo: image("ulgfsa.svg"),
  branch: "Engineering physics",
  academic-year: "2025-2026",
  tof: false,
  tot: false,
  toc: false,
  footer-text: "Homework"
)
#set par(justify: true)
#set page(paper: "a4")
#set text(size: 12pt)
#set page(numbering: "1")
#set math.equation(numbering: "(1)", number-align: bottom)
#import "@preview/zero:0.5.0": ztable
#import "@preview/theorion:0.4.1": *
#import cosmos.clouds: *
#show: show-theorion
#let theorem = theorem.with(fill: blue.lighten(85%))
#let theorem-box = theorem-box.with(fill: blue.lighten(85%))
#let theorem-box = theorem-box.with(radius: 5pt)
#import "@preview/physica:0.9.7": *
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Problem Statement

#quote-box[
#set text(fill: black)
Use the matched asymptotic method to provide a composite solution, accurate to
leading order plus one correction, of the following problem

$ epsilon y''(x) - 2(x^2+1) y'(x) + 1 = 0 $<eq:problem>

for $x in [0,1]$ with $y(0) = y(1) = 0$ and where $epsilon << 1$.

+ Develop asymptotic solutions. Create a composite solution.
+ For a given choice of parameters, compare your approximate solution with a
  numerical solution.
]

= Composite Solution

The @eq:problem is classified as singular because the small parameter $epsilon$ multiplies the highest‑order derivative. Consequently, setting $epsilon = 0$ reduces the order of the equation, making it impossible to satisfy both boundary conditions simultaneously. Hence, a rapid transition layer a boundary layer must develop near one of the endpoints.

== Location of the Boundary Layer

#theorem-box(title: "Boundary-layer location criterion")[
Consider the singularly perturbed linear differential equation: \
$ epsilon y''(x)+a(x) y'(x)+b(x) y(x)=f(x) $
with $x in [0,1]$ and the coefficient functions a,b,f are sufficiently smooth
Then the asymptotic structure of the solution exhibits a boundary layer at the endpoint(s) determined by the sign of ($a(x)$) as follows:
- If $a(x)>0 forall x in [0,1]$:,the layer is anchored at the left endpoint (x=0).
- If $a(x)<0 forall x in [0,1]$,the layer is anchored at the right endpoint (x=1).
- If $a(x)$ changes sign on $[0,1]$, a layer may develop in the interior of the domain, typically near the point where $a(x)=0$ or at the nearest endpoint where the characteristic direction reverses.
]
In @eq:problem, the coefficient of $y'$ is
$a(x) = -2(x^2 + 1)$.

Since $x^2 + 1 >= 1 > 0$ for all $x in [0,1]$, one has $a(x) < 0$ throughout
the domain. Therefore, the boundary layer is located near $x = 1$, and the
boundary condition $y(0) = 0$ is assigned to the outer region while $y(1) = 0$
is enforced by the inner solution.

== Outer Solution

Away from $x = 1$ all derivatives remain $O(1)$, consequently the term $epsilon y''$ is
negligible at leading order. The outer solution is sought as the two-term ansatz

$ y_"out"(x) = y_0(x) + epsilon y_1(x) + O(epsilon^2). $<eq:outer-ansatz>

Substituting @eq:outer-ansatz into @eq:problem and grouping terms by powers of
$epsilon$ yields separate problems at each order.

=== Leading order $O(epsilon^0)$

$ -2(x^2+1) y_0 '(x) + 1 = 0, quad y_0(0) = 0. $

Solving for $y_0 '$ and integrating:

$ y_0 '(x) = frac(1, 2(x^2+1)) quad ==> quad y_0(x) = frac(1,2) arctan(x) + A. $

Applying $y_0(0) = 0$ gives $A = 0$, yields

$ y_0(x) = frac(1,2) arctan(x). $<eq:y0>

=== First correction $O(epsilon^1)$

The term $epsilon y_0 ''$ carries over from the left-hand side of @eq:problem at
order $epsilon^1$, yielding

$ y_0 ''(x) - 2(x^2+1) y_1 '(x) = 0, quad y_1(0) = 0. $

Differentiating @eq:y0 twice:

$ y_0 '(x) = frac(1, 2(1+x^2)) quad ==> quad y_0 ''(x) = frac(-x, (1+x^2)^2). $

Substituting and solving for $y_1 '$:

$ y_1 '(x) = frac(-x, 2(1+x^2)^3). $

Integrating with the substitution $u = 1+x^2$, $dif u = 2x dif x$:

$ y_1(x) = frac(1, 8(1+x^2)^2) + B. $

Applying $y_1(0) = 0$ gives $B = -1/8$, producing

$ y_1(x) = frac(1, 8(1+x^2)^2) - frac(1,8). $<eq:y1>

The complete outer solution, valid away from $x = 1$, is therefore

$ y_"out"(x) = frac(1,2) arctan(x) + epsilon lr([frac(1, 8(1+x^2)^2) - frac(1,8)]) + O(epsilon^2). $<eq:outer>

#quote-box[
#set text(fill: black)
  *Remark.* The boundary condition $y(1) = 0$ is _not_ applied to the outer
  solution. It belongs to the inner (boundary layer) problem and will be used
  there.
]

== Boundary Layer Solution

=== Rescaling

Introduce the inner variable
$XX = frac(1-x, delta)$ such that $XX = O(1)$ inside the boundary layer.
$
y' = dv(y,x)= dv(YY(XX(x)),XX) dv(XX,x)= -1/delta YY '
$
where $YY(XX) equiv y(1 - epsilon XX)$ and primes for $YY$ now denote
derivatives with respect to $XX$.
In a analogous manner,one can deduce that

$
y'' =1/delta^2 YY ''
$

Substituting into @eq:problem:

$ epsilon/delta^2 YY '' + (4 - 4 epsilon XX + O(epsilon^2))1/delta YY ' +1 = 0, $<eq:inner-rescaled>

To have the first two terms of the same order we need
$ epsilon/delta^2 tilde 1/delta quad => quad delta tilde epsilon $
Finally,

$ YY '' + (4 - 4 epsilon XX + O(epsilon^2))YY ' + epsilon = 0, $<eq:inner-rescaledk>

By using a similar ansatz than the outer one,
$ YY(XX) = YY_0(XX) + epsilon YY_1(XX) + O(epsilon^2). $<eq:inner-ansatz>

Substituting @eq:inner-ansatz into @eq:inner-rescaledk and collecting by order:

==== Leading order $O(epsilon^0)$

$ YY_0 '' + 4 YY_0 ' = 0, quad YY_0(0) = 0. $

This is a first-order ODE in $YY_0 '$, with characteristic equation
$z^2 + 4z = 0$, giving roots $z = 0$ and $z = -4$. The general solution
satisfying $YY_0(0) = 0$ is

$ YY_0(XX) = C(e^(-4 XX) - 1), $<eq:y0-inner>

where the constant $C$ is determined by matching.

==== First correction $O(epsilon^1)$

Substituting @eq:inner-ansatz into @eq:inner-rescaledk at order $epsilon^1$, and
using $YY_0 ' = -4C e^(-4 XX)$:

$ YY_1 '' + 4 YY_1 ' = 4 XX YY_0 ' - 1 = -16 C XX e^(-4 XX) - 1, quad YY_1(0) = 0. $<eq:y1-inner>
#pagebreak()
A particular solution of the form
$YY_(1,p) = (D XX + E XX^2)e^(-4 XX) + F XX$
is substituted into @eq:y1-inner. Computing $YY_(1,p) '' + 4 YY_(1,p) '$
and matching coefficients:

#align(center)[
  #table(
    columns: 3,
    stroke: 0.5pt,
    [*Term*], [*Equation*], [*Result*],
    [$e^(-4 XX)$ (const)],  [$2E - 4D = 0$],   [$D = C$],
    [$XX e^(-4 XX)$], [$-8E = -16C$],    [$E = 2C$],
    [constant],                   [$4F = -1$],        [$F = -1 slash 4$],
  )
]

Adding the homogeneous solution $G e^(-4 XX) + H$ and applying
$YY_1(0) = 0$ gives $H = -G$, so

$ YY_1(XX) = (G + C XX + 2C XX^2) e^(-4 XX) - G - frac(XX,4). $<eq:y1-inner-sol>

The constant $G$ is again fixed by matching.

== Higher-Order Matching

The constants $C$ and $G$ are determined by requiring that the outer and inner
expansions agree to $O(epsilon)$ in the overlap region near $x = 1$.

=== Outer expansion near $x = 1$

Using Taylor series at $x = 1$ (writing $x - 1$ as small):

$ arctan(x) &= frac(pi,4) + frac(x-1,2) + O((x-1)^2), \
  frac(1,(1+x^2)^2) &= frac(1,4) - frac(x-1,2) + O((x-1)^2). $

Substituting into @eq:outer and retaining terms through $O(epsilon)$:

$ y_"out" approx underbrace(frac(pi,8), O(1)) + underbrace(frac(x-1,4), O(epsilon)) + epsilon lr([- frac(3,32) - frac(x-1,16)]). $<eq:outer-expanded>

=== Inner expansion as $XX -> +infinity$

As $XX -> +infinity$, all exponentials $e^(-4 XX) -> 0$. Using
$XX = (1-x)/epsilon$:

$ YY_"in" approx underbrace(-C, O(1)) + epsilon lr([- G - frac(XX,4)]) = -C - epsilon G + frac(x-1,4). $<eq:inner-expanded>

=== Determination of constants

Equating @eq:outer-expanded and @eq:inner-expanded at each order of $epsilon$:

#align(center)[
  #table(
    columns: 3,
    stroke: 0.5pt,
    [*Order*], [*Condition*], [*Result*],
    [$O(epsilon^0)$], [$-C = pi slash 8$],   [$C = -pi slash 8$],
    [$O(epsilon^1)$], [$-G = -3 slash 32$],  [$G = 3 slash 32$],
  )
]

The $O(epsilon)$ term proportional to $(x-1)$ matches automatically in both
expansions, confirming consistency.

The common asymptotic expansion in the overlap region is

$ y_"match" (x) = frac(pi,8) + frac(x-1,4) - frac(3 epsilon,32). $<eq:ymatch>

== Composite Solution

The uniformly valid composite solution is constructed as

$ y_"comp" (x) = y_"out" (x) + YY_"in" (XX) - y_"match" (x), $<eq:composite-formula>

where $y_"match"$ is subtracted once to avoid double-counting the common part
present in both $y_"out"$ and $YY_"in"$.

After substituting @eq:outer, @eq:y0-inner, @eq:y1-inner-sol, and @eq:ymatch
into @eq:composite-formula and simplifying, noting that
$e^(-4 XX) = e^(4(x-1) slash epsilon)$, the result is

$
y_"comp" (x) = frac(1,2) arctan(x)
  - frac(pi(2-x),8) e^(4(x-1) slash epsilon) \
  - frac(pi(1-x)^2, 4 epsilon) e^(4(x-1) slash epsilon)
  + epsilon lr([frac(1,8(1+x^2)^2) - frac(1,8) + frac(3,32) e^(4(x-1) slash epsilon)])
. $<eq:composite>

#quote-box[
#set text(fill: black)
  *Note on the term* $display(frac(pi(1-x)^2, 4 epsilon))$.
  Although this factor appears to grow as $epsilon -> 0$, inside the boundary
  layer $(1-x) = O(epsilon)$, so $(1-x)^2 / epsilon = O(epsilon)$. The term
  is therefore genuinely $O(epsilon)$ in the BL and transcendentally small
  (exponentially suppressed) everywhere else.
]

=== Verification of boundary conditions

- *At $x = 0$:* one has $e^(4(0-1) slash epsilon) = e^(-4 slash epsilon)$, which
is transcendentally small for any fixed $epsilon > 0$. Consequently all
exponential terms vanish and

$ y_"comp" (0) approx frac(1,2) arctan(0) + epsilon lr([frac(1,8) - frac(1,8)]) = 0. $

- *At $x = 1$:* one has $e^0 = 1$ and $(1-x)^2 = 0$, so

$ y_"comp" (1) = frac(pi,8) - frac(pi,8) - 0 + epsilon lr([frac(1,32) - frac(1,8) + frac(3,32)]) = 0 + epsilon dot 0 = 0. $

Both boundary conditions are satisfied exactly.

= Numerical Results

The composite solution @eq:composite is compared against a numerical reference
obtained by solving the boundary value problem with `scipy.integrate.solve_bvp()`
at tolerance $10^(-10)$. The second-order ODE is rewritten as the first-order
system $y_1' = y_2$, $y_2' = (2(x^2+1) y_2 - 1)/epsilon$, with boundary
conditions $y_1(0) = 0$ and $y_1(1) = 0$.The $epsilon$ is considered at 0.2 and 0.01.
#figure(
  image("figures/plot_0.2.svg", width: 100%),
  caption: [
    Comparison of the numerical reference (blue), the $O(epsilon)$ composite
    solution (orange), the outer-only solution (green dotted), and the
    inner-only solution (red dash-dotted), for $epsilon= 0.2$.
  ]
)<fig:comparison>

@fig:comparison and @fig:comparison2 illustrates the behaviour of each asymptotic approximation across the domain.

- *Outer solution.* The outer approximation captures the slow variation of $y$
  across most of $[0,1]$, but diverges from the true solution in the vicinity of
  $x = 1$, where it is unable to satisfy the right boundary condition.

- *Inner solution.* The inner solution resolves the rapid variation near $x = 1$
  and satisfies $y(1) = 0$ by construction. Away from the boundary layer it
  approaches a non-zero constant (its matching limit $-C = pi/8$), so it does
  not constitute a valid approximation over the full domain on its own.

- *Composite solution.* The composite correctly captures both the smooth outer
  behaviour and the sharp boundary layer.Due to $epsilon = 0.2$ a small discrepancy appears near the maximum of the curve, which is expected since the asymptotic ordering
  $epsilon << 1$ becomes less rigorous.

#figure(
  image("figures/plot_0.01.svg", width: 100%),
  caption: [
    Comparison of the numerical reference (blue), the $O(epsilon)$ composite
    solution (orange), the outer-only solution (green dotted), and the
    inner-only solution (red dash-dotted), for $epsilon= 0.01$.
  ]
)<fig:comparison2>
