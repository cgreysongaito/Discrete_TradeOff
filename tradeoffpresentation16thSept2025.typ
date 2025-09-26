#import "@preview/thmbox:0.3.0": *

#show: thmbox-init(counter-level: 1)


#import "@preview/slydst:0.1.4": *
#import "@preview/wrap-it:0.1.1": wrap-content

#show: slides.with(
  title: "Adding a fecundity-survival trade-off to a discrete population model with maturation delay",
  subtitle: "Chris Greyson-Gaito, Sabrina Streipert, Gail Wolkowicz",
  date: "September 18th, 2025",
  layout: "medium",
  ratio: 16/9,
  title-color: none,
)
#set text(16pt)

== Delay Recurrence Population Models
#v(1cm)
#align(center)[
$X_(t+1)=F(X_t,X_(t-tau))$

blends

simplest scalar form $X_(t+1)=F(X_t)$

with

age-structure form $arrow(X)=arrow(F)(arrow(X)_t)$
]

== Problem?
#v(2cm)
Many delay recurrence models have not accounted for the losses that occur between birth and maturity.
#v(1cm)
*Streipert & Wolkowicz 2023 accounted for losses during maturation delay.*



== Streipert & Wolkowicz 2023 (1)
Derivation:

$X_(t+1)=p(X_t)X_t+p(X_t)g(X_t)X_t= F_1(X_t)+F_2(X_t)$

where 
- $F_1$ is the number of inviduals in the population alive and sexually mature at the end of the interval from $t-1$ to $t$ that survive the time interval from $t$ to $t+1$
- $F_2$ describes the density of individuals that were born at the beginning of the breeding cycle from $t$ to $t + 1$ that reach maturity at the end of this cycle (i.e. at time t + 1)

== Streipert & Wolkowicz 2023 (2)
#v(1cm)
We could make the survival probability different between mature and immature individuals.

$X_(t+1)=p(X_t)X_t+hat(p)(X_t)g(X_t)X_t$

This still assumes that newborn individuals join the population and procreate at the beginning of the breeding cycle immediately following the one in which they were born.


== Streipert & Wolkowicz 2023 (3)
#v(1cm)
Many species reproduce after multiple breeding cycles

$X_(t+1)=p(X_t)X_t+g(X_(t-tau))tilde(p)(tau,X_(t-tau))X_(t-tau)$

where

- $tau+1$ is the number of breeding seasons from birth until they reach sexual maturity
- we assume that immature individuals interact with only other immature individuals of the same cohort
== Streipert & Wolkowicz 2023 (4)
#v(0.5cm)
*Example - Beverton-Holt*

Without delay:

$X_(t+1) = (1+g(t))/(1+q(t))X_t = (1+r)/(1+d+c X_t)X_t$


$X_(t+1)=(1+r)/(1+d+c X_t)X_t = 1/(1+d+c X_t)X_t + 1/(1+d+c X_t)r X_t$

Now with a delay:

$X_(t+1)=(1+r)/(1+d+c X_t)X_t = 1/(1+d+c X_t)X_t + tilde(p)(tau,X_(t-tau))r X_(t-tau)$

== Streipert & Wolkowicz 2023 (5)

*Example - Beverton-Holt*

$X_(t+1)=(1+r)/(1+d+c X_t)X_t = 1/(1+d+c X_t)X_t + tilde(p)(tau,X_(t-tau))r X_(t-tau)$

#v(1cm)
What is $tilde(p)(tau,X_(t-tau))$?
#v(1cm)

Through solving via recurrence:

$tilde(p)(tau,X_(t-tau))= D/ (D(1+D)^(tau+1)+((1+D)^(tau+1)-1)C r X_(t-tau))$




== Streipert & Wolkowicz 2023 (6)
#v(2cm)
#align(center)[
#text(18pt)[*Increasing the maturation delay decreased the positive interior equilibrium until extinction occurred.*]
]
== But, maturation delays are very common in species.

#align(center)[#image("maturationdelay.png", width: 50%)]

#align(center)[Budd et al. 2024]

== Life history trade-offs!
Negative association between two traits

e.g. reproduction, survival, reserves ...

#align(center)[#image("energytradeoff.png", width: 60%)]

#align(center)[#text(12pt)[Zera and Harshman 2001]]

== Our trade-off
#align(center)[
  #v(2cm)
*Longer maturation delay*

Lower total survival

Higher fecundity

]

== Fecundity

#align(center)[#image("fecundityimage.png")]

== Our  models

#align(center)[#image("figs/schematic_draft1.png", width: 68%)]

== Scenario i) Constant survival of immature individuals
#v(3cm)
1) BevertonHolt  - Constant

2) Ricker - Constant

== BevertonHolt - Constant
#align(center)[
#image("figs/BevHoltI_equi.png", width: 60%) 

$overline(p)$ is the constant survival fraction of immature individuals

]

== Ricker - Constant (even $tau$)
#v(1.5cm)
#align(center)[#image("figs/RickerConstanttaueven.svg")]

== Ricker - Constant ($tau=3$)
#v(1.5cm)
#align(center)[#image("figs/RickerConstanttau3.png")]

== Ricker - Constant (two par)
#v(1.5cm)
#align(center)[#image("figs/RickerConstant_apbifurcation_wtauorbits.svg")]

== Scenario ii) Cohort-density dependent survival of immature individuals
#v(2cm)
1) BH - BH

2) Ricker - Ricker

== BH - BH

#align(center)[#image("figs/BevHoltBevHolt_tauequi.svg")]

== Ricker - Ricker (1)
#v(1cm)
#h(3cm)#text(14pt)[Immature density dependence#h(1.5cm)Mature density dependennce]
#align(center)[#image("figs/aCbetabifurcation_RickerRicker_wtauorbitsa.png", width:100%)]

== Ricker - Ricker (2)
#align(center)[#image("figs/aCbetabifurcation_RickerRicker_wtauorbitsb.png")]

== Conclusion

- Adding the fecundity-survival trade-off lead to an optimal maturation delay
-- Consistent with the biological literature
#v(0.8cm)
- Fecundity-survival trade-off affected population stability (with Ricker survival function)
-- Increasing maturation delay:

No interior equilibrium $->$ Stable interior equilibrium $->$ 

#h(2cm)Oscillations $->$ Remove oscillations $->$ Extinction

== Further extensions
#v(2cm)
#text(22pt)[
- Use distributed delay instead of discrete delay
- Scale loss rates with maturation delay
- Explore other life history trade-offs
]

== Questions?

#align(center)[#text(20pt)[*Thank you!*]


#box(height: 25pt, image("logo-purple.svg", width: 5%))cgreysongaito\@ecoevo.social

#v(1cm)

#box(height: 20pt, image("64px-Bluesky_Logo.svg.png", width: 5%))\@cgreysongaito.bsky.social

#v(1cm)

Made in #box(height: 14pt, image("typst.svg"))

]

// #wrap-content(text[Made in],figure(image("typst.svg")))
// ]


