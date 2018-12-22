Statistical Rethinking book Errata
==========

page 42: Just below R code box 2.6, the text says that map requires a list of start values. It does not, as long as priors are provided for each parameter. And indeed the example in box 2.6 does not contain a list of start values.

page 83: Overthinking: Model definition to Bayes' theorem again: first term of the numerator and denominator might best be parenthesized to unambiguously discourage factoring in the prior once for each observation (i times), i.e.,
<a href="https://www.codecogs.com/eqnedit.php?latex=${\left(\Pi_i\text{Normal}(h_i\mid\mu,\sigma))&space;\text{Normal}(\mu\mid&space;178,20)&space;\text{Uniform}(\sigma\mid&space;0,50)\right)&space;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?${\left(\Pi_i\text{Normal}(h_i\mid\mu,\sigma))&space;\text{Normal}(\mu\mid&space;178,20)&space;\text{Uniform}(\sigma\mid&space;0,50)\right)&space;}" title="${\left(\Pi_i\text{Normal}(h_i\mid\mu,\sigma)) \text{Normal}(\mu\mid 178,20) \text{Uniform}(\sigma\mid 0,50) }" /></a>

page 87: The marginal description of the model reads "mu ~ dnorm(156, 10)" but the model is Normal(178, 20). Same error on p 95 and in code 4.38. It is corrected in code 4.39.

page 95-96: dnorm(156,100) should be dnorm(178,100) in both model presentation and then R code on top of page 96.

page 103, R code 4.50: The ``post`` object implied here is the one from R code 4.46: ``post <- extract.samples(m4.3)``.

page 125: Below R code 5.4, "The posterior mean for age at marriage, ba, ..." 'ba' should be 'bA'.

page 156, near top: "In fact, if you try to include a dummy variable for apes, you'll up with..." Should be "you'll end up with".

page 196-200: The data.frame d has 17 cases. However in the discussion of the four models (on e.g. page 200), the text repeatedly refers to 12 cases.

page 212, the next-to-last sentence on the page refers to "the Rethinking box at the end of this section." That box is not in the text.

page 237 Exercise H1: "...index variable, as explained in Chapter 6.",
should be chapter 5 (at least that's their first appearance)

page 253 ("...the functions postcheck, link and sim work on map2stan
just as they do on map models...") postcheck appears somewhat out of thin air. Need a better introduction to it.

page 314: "Islands that are better network acquire or sustain more tool types.": network should be networked.

page 331, line 1: "a women" -> "a woman"

page 435: "FIGURE 14.4 display ... imputed neocortex values in blue ...
shown by the blue line segments". The imputed values are actually the
open black dots (and corresponding black line segments) as the caption
of the figure correctly states.

