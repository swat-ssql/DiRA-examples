**Draft version, not for distribution. Please do not share this link.**

The core idea behind the DiRA package is that once you can visualize and interact with a regression surface in 3D, you can productively examine and interpret directions within the regression surface. The DiRA package specifically allows you to:

  - Visualize a regression surface and confidence intervals for a model of the form y ~ x1 + x2 or y ~ x1 + x2 + x1 * x2.
  - Visualize the scatterplot of observations measured by y, x1, and x2 along with the regression surface.
  - Visualize regression surfaces estimated from a categorical explanatory variable with three categories.
  - Visualize specifics directions in a regression surface, including
  
      * The marginal effects and confidence intervals of x1 and x2 in the regression surface.
      * The direction defined by an omitted variable model, which is a type of nested model, y ~ x1 or y ~ x2. This direction can show the confidence intervals associated with the nested model.
      * Any rotations of the original model. This is useful for examining the confidence intervals associated with the regression surface for the direction defined by a nested model. It is also useful for examining all possible pairwise comparisons for a 3 level categorical variable, examining intersectionality from 2 categorical variables, or for examining principal component directions.
      * Concisely view the estimated effect and standard errors of the original regression model, as well as any nested or rotated versions of the original regression model.


This repo was initially generated from a bookdown template available here: https://github.com/jtr13/bookdown-template.

Demo Video
A demo video showing how to create a bookdown book following these instructions: http://bit.ly/fiveminutebookdown

Additional features
Please consult the official guide to bookdown: https://bookdown.org/yihui/bookdown
