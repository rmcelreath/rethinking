# Statistical Rethinking book Errata

## 2nd Edition

page 41: Should be Figure 2.6

page 108 (R code 4.61): mu.HPDI is not defined a this point - I presume you mean mu.PI (defined in R code 4.56 on page 105)?

page 121 (Problem 4M6): Misleading to state that the variance is 64 cm. If any units are used, it should be $cm^2$

page 127 (Fig 5.3 text): implausibly

page 130: Last paragraph: "And D and M are associated with another, because A influences them both"

page 132: Fourth paragraph (after the formula): "but here I've chosen M for marriage rate"

page 148: First paragraph: "But if have..." -> "But I have"

page 151-152: R code 5.41 actually produces the lower-right plot in Figure 5.9

page 213: Overthinking box: function name is sim_train_test now in 2nd ed.

page 228: R code 7.28: Should be 40.9, to tie the number more tightly to the paragraph above

page 248: R code 8.14: Should be precis(m8.3, depth=2 )

page 273: Last equation, should end with "log p(\mu_x | 0, 0.5)" (substitute vertical bar for the comma)

page 278: R code 9.10: Beware that the real HMC2 function also defines dH in the returned list. This is used in the plot (R code 9.7). It is the difference in total energy (proposed_U + proposed_K) - (current_U + current_K).

page 331: last paragraph: should be CONTRASTS.

page 338: "results for the first two chimpanzees" should really be (IMHO): "partial results for the first two treatments (the total data frame has 7x4=28 rows, one each per actor and treatment combination."

page 339: First equation: typo in the numerator - should be 9! (in accordance with the binomial coefficient formula)

page 343 (R code 11.32): dat_list$dept_id <- coerce_index(d$dept)   # use coerce_index and the original data to construct new index variable

page 378: First line, should be "probability that the monks did drink"

page 410: Mismatch between R.code 13.7 where a_bar is 1.5 and the second paragraph, where it is assumed to be 1.4

page 453: (Overthinking): Typo in paragraph 4: Should be compose_noncentered.

page 463: Typo in para 5, should be: "when one household gives more, the other gives less"

page 469: Typo in para 3, should be: "famililar Poisson probability"

page 472: Type in para 1, should be: "those k parameters"

page 475: Typo in Rcode 14.45, should be: exp( post$a + post$b*lp )

page 477: Typo, first line in ch.14.5.2, should be "more or less distant"

## 1st Edition

page 13: "What does mean to take a limit..." is missing the word "it".

page 42: Just below R code box 2.6, the text says that map requires a list of start values. It does not, as long as priors are provided for each parameter. And indeed the example in box 2.6 does not contain a list of start values.

page 66, end of first paragraph: 'the right-hand plot' should be 'the bottom plot'.

page 76, Overthinking box, first paragraph: "You're computer already knows it" should read "Your computer..."

page 87: The marginal description of the model reads "mu ~ dnorm(156, 10)" but the model is Normal(178, 20). Same error on p 95 and in code 4.38. It is corrected in code 4.39.

page 95-96: dnorm(156,100) should be dnorm(178,100) in both model presentation and then R code on top of page 96.

page 103, R code 4.50: The ``post`` object implied here is the one from R code 4.46: ``post <- extract.samples(m4.3)``.

page 125: Below R code 5.4, "The posterior mean for age at marriage, ba, ..." 'ba' should be 'bA'.

page 156, near top: "In fact, if you try to include a dummy variable for apes, you'll up with..." Should be "you'll end up with".

page 196-200: The data.frame d has 17 cases. However in the discussion of the four models (on e.g. page 200), the text repeatedly refers to 12 cases.

page 212, the next-to-last sentence on the page refers to "the Rethinking box at the end of this section." That box is not in the text.

page 215, first paragraph: "despite it's plausible superiority" should be "despite its plausible superiority".

page 237 Exercise H1: "...index variable, as explained in Chapter 6.",
should be chapter 5 (at least that's their first appearance)

page 253 ("...the functions postcheck, link and sim work on map2stan
just as they do on map models...") postcheck appears somewhat out of thin air. Need a better introduction to it.

page 314: "Islands that are better network acquire or sustain more tool types.": network should be networked.

page 331, line 1: "a women" -> "a woman"

page 386, problem 12H1, first paragraph: 'By the year 200' should read 'By the year 2000'.

page 403: The average effect in the P *C interaction model is typed βP but I think should be βPC.

page 435: "FIGURE 14.4 display ... imputed neocortex values in blue ...
shown by the blue line segments". The imputed values are actually the
open black dots (and corresponding black line segments) as the caption
of the figure correctly states.

