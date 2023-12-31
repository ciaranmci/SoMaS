History Matching
================
Ciarán McInerney
2023-06-27

``` r
# Load packages.
library(tidyverse)
library(FME)
library(RobustGaSP)
library(gridExtra)
```

## TL;DR

**History-matching methods are just fancy ways to get rid of input
values that give ridiculous output values. You do this by considering
what kind of previously-used input values gave ridiculous output values,
and deciding not to carry them forward.** <br/><br/> **If you’ve ever
made a meal and added too much salt, you will have learned how much salt
makes for a bad meal. Next time you make the meal and you are unsure
about how much salt to add, you will consider your previous (a.k.a.
historical) salting and say “*Well, there is definitely no point in
adding more salt than I did last time*”. Congratulations, you have just
applied history matching to constrain the range of inputs (i.e. salt) so
that the outputs (i.e. tastiness of your meal) are more likely to be
tasty / desirable / sensible / realistic!** <br/><br/> **We apply
history matching to constrain simulation and emulation models because
computer models will do whatever calculations you program them to, even
if the program, the inputs, and the outputs are ridiculous. Never
forget: Computers aren’t smart, they do dumb things fast. Because
simulators are models of reality, we constrain them using theory and
historical observations of reality. Because emulators are models of
models, we constrain them using historical observations of the
underlying model that is being emulated.** <br/><br/> *(Side note:
Incidentally, check out this [Art of Manliness podcast episode with
Daniel
Holzman](https://www.artofmanliness.com/living/food-drink/podcast-812-chef-vetted-answers-to-your-cooking-faqs/)
to learn a lot about appropriate salting and more.)* <br/><br/>
<br/><br/>

## Why?

So, what is the motivation behind history matching? <br/><br/> Well, no
body likes getting things wrong. But because we can’t know everything,
we get through life by developing heuristics and models of the way the
world works. These heuristics and models help us to choose the actions
we take that lead to the outcomes we perceive. Of course, the future is
emergent and uncertain so we can’t guarantee particular outcomes.
Instead, we must settle for uncertainty. Or must we? <br/><br/> Well,
yes, we do, but could we *improve* our certainty about outcomes? In
other words, could we improve what we know about the outputs for our
given inputs? The options that are open to us are either:

1.  maximise the fidelity of our model of how the world works, on the
    assumption that its our messy model that results in messy outputs,
    or.

2.  narrow our range of inputs, on the assumption that the less variety
    going in results in less variety coming out.

Option 2 is the essence of history matching and seems very sensible.
After all, what’s the point in doing things that don’t result in the
outcomes you desire? As someone wise once said, [doing the same thing
over and over and expecting different results is the definition of
insanity](https://quoteinvestigator.com/2017/03/23/same/). <br/><br/>
<br/><br/>

## Seriously, what?

You can think of history matching as optimising performance over many
iterations, but only by changing what you put in to a system / model
rather than trying to change the system / model. It is just like
learning to win at a game: the rules are set (i.e. the model is fixed)
but while playing the game over and over again, you will try a variety
of tactics and eventually learn what leads to success and what leads to
failure. The keys are to stay humble, determined, and not to stubbornly
stick to worn-out ideas that haven’t born fruit. <br/><br/> To bridge
the gap between this prosaic explanation and a mathematically-grounded
understanding, I found [Williamson et
al.’s](https://www.sci-hub.wf/10.1007/s00382-013-1896-4) idea of the Not
Ruled Out Yet space to be useful. The Not Ruled Out Yet space is the
space of inputs (i.e. the scope of possible values for the inputs) that
remains after I’ve gotten rid of the input values that produce
ridiculous and unlikely outputs. We can come up with any number of ways
to decide what should be ruled in and ruled out. [Adrianakis et
al.](https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1003968&type=printable)
have a similar idea called the “*non-implausible region*”. <br/><br/>
For example, the implausibility measure used by [Williamson et
al.](https://www.sci-hub.wf/10.1007/s00382-013-1896-4) is the difference
between average predicted output value and the observed output values,
scaled to the standard deviation of all these differences:

$$
I(\ parameter\ value_{j}\ ) = max\{\ I_{i}(parameter\ value_{j})\ \}
$$

where

$$
I_{i}(parameter\ value_{j}) = \frac{ |\ observation_{i} - E[\ prediction_{i}( parameter\ value_{j})\ ]\ | } {\sqrt{ Var[\ observation_{i} - E[\ prediction_{i}(parameter\ value_{j})\ ]\ ] }}
$$

and where $i$ refers to runs of the model, and $j$ refers to inputs for
particular parameter values. You use the implausability measure to rank
parameter values then trim off the worst. <br/><br/> Any choice of rule
comes with trade-offs. The implausibility measure, for example, can be
compared across different predictive models even with different outputs
because it is a unitless, scaled score. But, if the process of scoring,
ranking, and triming always trims the most-extreme scores, then repeated
application (called “refocussing”) can result in a situation like
[Zeno’s dichotomy
paradox](https://en.wikipedia.org/wiki/The_Indefatigable_Frog) and can
result in the Not Ruled Out Yet space vanishig. To prevent this, some
absolute, minimum size of the Not Ruled Out Yet space needs to be
set…but couldn’t we just have set an absolute, maximum acceptable
difference in the first place? <br/><br/> As an alternative perspective
to [Williamson et
al.’s](https://www.sci-hub.wf/10.1007/s00382-013-1896-4) Not Ruled Out
Yet space, you can think of history matching as Monte Carlo filtering.
This is what sensitivity-analysis guru Andrea Saltelli summarise in
[Saltelli et al.](https://sci-hub.wf/10.1002/9780470725184) (2008, page
248):

> “one samples the space of the input factors as in ordinary \[Monte
> Carlo\], and then categorizes the corresponding model output as either
> within or without the target region (the terms behaviour, $B$, or
> nonbehaviour, $\bar{B}$, are used). This categorization is then mapped
> back onto the input factors, each of which is thus also partitioned
> into a behavioural and nonbehavioural subset.”

<br/><br/>

## How?

I’ve already explained how history matching is sometimes done -
e.g. scoring with the implausibility measure, ranking, then trimming -
but there are many ways. Earlier, I said history matching was like
optimising performance, so you might not be surprised to hear that there
are many history-matching methods that use optimisation routines on data
from feedback systems. Another approach is to start with some idea of
what sensible inputs are and then update your understanding as time goes
on, i.e. the Bayesian approach. I’ll briefly present some of these,
below. First, a quick demo. <br/><br/>

### Matching on limits.

Let’s take a simple and deterministic model that accepts $x$ as a single
input and outputs $y=f(x)$. Let’s also pretend that history has taught
us that outputs beyond $y=f(x)=0\pm0.5$ are complete nonsense. The model
might produce these output values, but they are unacceptable to us. I’ll
use Latin Hypercube sampling to provide a good spread of input values,
and we’ll figure out which ones lead to acceptable outputs and which
ones don’t. <br/><br/> Before we do any calculations to determine the
unacceptability of our inputs, we can (in our simple example) check a
simple plot.

``` r
# Define a simple model.
mod <- function(x)
{
  df <- data.frame(input = x,
                   output = 2*x + 3*x*sin(5*pi*(x-0.1)/0.4)
                   ) %>% dplyr::arrange()
  colnames(df) <- c("input", "output")
  return(df)
}
# Set input values.
possVals <- data.frame(min = 0, max = 1)
rownames(possVals) <- c("input")
n_obs <- 50
df_input <- FME::Latinhyper(possVals, n_obs) %>% data.frame() %>% dplyr::arrange()
# Calculate the output.
df_output <- mod(df_input)
# Visual check of acceptability.
(p_simple_HM <-
    df_output %>%
    ggplot() +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = -0.5, ymax = 0.5, fill = "palegreen", alpha = 0.2) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = Inf, fill = "lightcoral", alpha = 0.2) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -0.5, fill = "lightcoral", alpha = 0.2) +
    geom_hline(yintercept = 0.5) + geom_hline(yintercept = -0.5) +
    geom_line(aes(x = input, y = output)) +
    geom_point(aes(x = input, y = output), size = 5, shape = 1) +
    scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 20))
)
```

<img src="HistoryMatching_explainer_files/figure-gfm/Simple history matching plot-1.png" width="50%" style="display: block; margin: auto;" />
Oh dear! Many points are in our red unacceptable region, across almost
the entire range of inputs. Let’s now name-and-shame the inputs that are
leading to unacceptable outputs.

``` r
# Define our rule for unacceptable outputs.
calculate_unacceptability <- function(y)
{
  unacceptable <- y %>% dplyr::filter(output > 0.5 | output < -0.5) %>% dplyr::arrange(input)
  acceptable <- y %>% dplyr::anti_join(unacceptable) %>% dplyr::arrange(input)
  acceptability_check <- list(unacceptable, acceptable)
  names(acceptability_check) <- c("unacceptable", "acceptable")
  
  return(acceptability_check)
}
# Get the acceptable and unacceptable inputs.
yes_and_no <- calculate_unacceptability(df_output)
# Update the previous plot.
(p_HM_shame <-
  ggplot(yes_and_no$unacceptable) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -0.5, ymax = 0.5, fill = "palegreen", alpha = 0.2) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = Inf, fill = "lightcoral", alpha = 0.2) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -0.5, fill = "lightcoral", alpha = 0.2) +
  geom_hline(yintercept = 0.5) + geom_hline(yintercept = -0.5) +
  geom_point(data = yes_and_no$unacceptable, aes(x = input, y = output), color = "red4", shape = 4) +
  geom_point(data = yes_and_no$acceptable, aes(x = input, y = output), size = 5, shape = 1) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 20))
)
```

<img src="HistoryMatching_explainer_files/figure-gfm/Simple history matching calculation-1.png" width="50%" style="display: block; margin: auto;" />
Unfortunately for our “simple” function, the history-matching process
leaves it to us to handle discontinuous blocks of acceptable inputs. But
this is a contrived example. Often, we don’t know the underlying
data-generating mechanism and are tasked with fitting a statistical
model to our observations. In such a context, history matching really
shines because we use it to filter the inputs before refitting our
model. <br/><br/> Below, I fit a Gaussian process emulator to our
observations before and after using history matching to filter for
acceptable inputs. (You can find out more about Gaussian process
emulators in [my other blog](https://github.com/ciaranmci/SoMaS) )

``` r
# Train emulator on all observed data.
GPE_model_preHM <- RobustGaSP::rgasp(design = df_output$input,
                                     response = df_output$output)
# Get emulator predictions.
input_plot_vals <- data.frame(inputs = seq(0, 1, 0.01))
df_GPE_model_preHM <-
  predict(GPE_model_preHM, input_plot_vals) %>%
  as.data.frame() %>%
  dplyr::bind_cols(input_plot_vals)
# Make first plot of Gaussian process emulated outputs.
p_GPE_preHM <-
  df_GPE_model_preHM %>%
  ggplot() +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -0.5, ymax = 0.5, fill = "palegreen", alpha = 0.2) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = Inf, fill = "lightcoral", alpha = 0.2) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -0.5, fill = "lightcoral", alpha = 0.2) +
  geom_ribbon(aes(x = inputs, ymin = lower95, ymax = upper95), fill = "grey70") +
  geom_point(data = df_output, aes(x = input, y = output), size = 5, shape = 1) +
  geom_line(aes(x = inputs, y = mean), linewidth = 1) +
  labs(title = 'Emulator fitted to all values',
       subtitle = 'Note the practically-invisible uncertainty in the\nemulators mean prediction (black line)\nwhen all observations are used.') +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5), limits = c(-1,5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 10))

# Train emulator on all observed data.
GPE_model_postHM <- RobustGaSP::rgasp(design = yes_and_no$acceptable$input,
                                      response = yes_and_no$acceptable$output)
# Get emulator predictions.
df_GPE_model_postHM <-
  predict(GPE_model_postHM, input_plot_vals) %>%
  as.data.frame() %>%
  dplyr::bind_cols(input_plot_vals)
# Make second plot of Gaussian process emulated outputs.
p_GPE_postHM <-
  df_GPE_model_postHM %>%
  ggplot() +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -0.5, ymax = 0.5, fill = "palegreen", alpha = 0.2) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = Inf, fill = "lightcoral", alpha = 0.2) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -0.5, fill = "lightcoral", alpha = 0.2) +
  geom_ribbon(aes(x = inputs, ymin = lower95, ymax = upper95), fill = "grey70") +
  geom_point(data = yes_and_no$acceptable, aes(x = input, y = output), size = 5, shape = 1, color = "green4") +
  geom_point(data = yes_and_no$unacceptable, aes(x = input, y = output), size = 5, shape = 4, color = "red4") +
  geom_line(aes(x = inputs, y = mean), linewidth = 1) +
  labs(title = 'Emulator fitted to acceptable values',
       subtitle = 'Note how the emulator\'s uncertainty increased\nin the area with the many unacceptable\nvalues.') +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5), limits = c(-1,5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 10))
# Final plot.
gridExtra::grid.arrange(p_GPE_preHM, p_GPE_postHM, ncol = 2)
```

<img src="HistoryMatching_explainer_files/figure-gfm/GPE before and after HM-1.png" width="80%" style="display: block; margin: auto;" />
As I note in the subplot headers, this particular emulator is uncertain
about outputs in the regions where history matching removed unacceptable
inputs. Nevertheless, the mean emulation (the black line) has been
pulled closer to the acceptable range of inputs, which is what we wished
for. Unfortunately, we are left with a lo tof uncertainty. We should be
careful what we wish for. <br/><br/> <br/><br/>

### Relative measures of uncertainty.

In the previous example, we optimised our model by filtering out inputs
that led to unacceptable outputs. What we meant by “unacceptable” was
based on absolute values of the output. Alternatively, we might a
relative measure of uncertainty. <br/><br/> Let’s use the same function
as in the previous example. We have already fit an emulator model to all
data, so we will use that for our predictions (or emulations, to be more
correct). To “collect” observations, we will calculate the outputs for
the whole range of inputs using our deterministic function. We will plot
the acceptable and unacceptable regions as we iteratively trim the
most-unacceptable 80% using the implausibility measure that I described
earlier. This process is called refocussing.

``` r
# Set data.
trim_percentile <- 0.2
df_refocus <- df_output
observed_output <- mod(input_plot_vals)$output
# Make plot of all data in "acceptable" zone, with line for the emulated output.
p_refocus_1 <-
    df_GPE_model_preHM %>%
    ggplot() +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "palegreen", alpha = 0.2) +
    geom_ribbon(aes(x = inputs, ymin = lower95, ymax = upper95), fill = "grey70") +
    geom_line(aes(x = inputs, y = mean)) +
    geom_point(data = df_output, aes(x = input, y = output), size = 5, shape = 1) +
    scale_y_continuous(breaks = seq(-1,5,1), limits = c(-2, 5)) + labs(y = "outputs") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 10))
# Calculate the implausibility for each input.
diff <- observed_output - df_GPE_model_preHM$mean
imp = abs(diff) / sd(diff)
trim_pt <- quantile(imp, probs = trim_percentile)[[1]]
# Define the boundaries of the Not Rule Out Yet (NROY) region.
NROY_region <- dplyr::bind_cols(observed_output = observed_output, imp = imp) %>% dplyr::filter(imp < trim_pt) %>% dplyr::select(observed_output) %>% dplyr::summarise(y_max = max(observed_output), y_min = min(observed_output))
print(NROY_region)
# Filter the training data to observations within the NROY region.
trim_emul_input <- df_refocus %>% dplyr::filter(output <= NROY_region$y_max & output >= NROY_region$y_min)
df_refocus <- df_refocus %>% dplyr::mutate(trim1 = ifelse(output <= NROY_region$y_max & output >= NROY_region$y_min, 0, 1))
# Train new emulator on acceptable inputs only, then predict.
updated_GPE <- RobustGaSP::rgasp(design = trim_emul_input$input,
                                response = trim_emul_input$output)
df_updated_GPE_preds <-
  predict(updated_GPE, input_plot_vals) %>%
  as.data.frame() %>%
  dplyr::bind_cols(input_plot_vals)
# Make first trimmed plot.
p_refocus_2 <-
    df_updated_GPE_preds %>%
    ggplot() +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = NROY_region$y_min, ymax = NROY_region$y_max, fill = "palegreen", alpha = 0.2) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = NROY_region$y_max, ymax = Inf, fill = "lightcoral", alpha = 0.2) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = NROY_region$y_min, fill = "lightcoral", alpha = 0.2) +
    geom_ribbon(aes(x = inputs, ymin = lower95, ymax = upper95), fill = "grey70") +
    geom_line(aes(x = inputs, y = mean)) +
    geom_point(data = df_refocus %>% dplyr::filter(trim1 == 0), aes(x = input, y = output), size = 5, shape = 1) +
    geom_point(data = df_refocus %>% dplyr::filter(trim1 == 1), aes(x = input, y = output), size = 5, shape = 4, color = "red4") +
    scale_y_continuous(breaks = seq(-1,5,1), limits = c(-2, 5)) + labs(y = "outputs") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 10))
# Calculate the implausibility for each input, again.
diff <- observed_output - df_updated_GPE_preds$mean
imp = abs(diff) / sd(diff)
trim_pt <- quantile(imp, probs = trim_percentile)[[1]]
# Define the boundaries of the Not Rule Out Yet (NROY) region.
NROY_region <- dplyr::bind_cols(observed_output = observed_output, imp = imp) %>% dplyr::filter(imp < trim_pt) %>% dplyr::select(observed_output) %>% dplyr::summarise(y_max = max(observed_output), y_min = min(observed_output))
print(NROY_region)
# Filter the training data to observations within the NROY region.
trim_emul_input <- df_refocus %>% dplyr::filter(output <= NROY_region$y_max & output >= NROY_region$y_min)
df_refocus <- df_refocus %>% dplyr::mutate(trim2 = ifelse(output <= NROY_region$y_max & output >= NROY_region$y_min, 0, 1))
# Train new emulator on acceptable inputs only, then predict.
updated_GPE <- RobustGaSP::rgasp(design = trim_emul_input$input,
                                response = trim_emul_input$output)
df_updated_GPE_preds <-
  predict(updated_GPE, input_plot_vals) %>%
  as.data.frame() %>%
  dplyr::bind_cols(input_plot_vals)
# Make second trimmed plot.
p_refocus_3 <-
    df_updated_GPE_preds %>%
    ggplot() +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = NROY_region$y_min, ymax = NROY_region$y_max, fill = "palegreen", alpha = 0.2) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = NROY_region$y_max, ymax = Inf, fill = "lightcoral", alpha = 0.2) +
    geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = NROY_region$y_min, fill = "lightcoral", alpha = 0.2) +
    geom_ribbon(aes(x = inputs, ymin = lower95, ymax = upper95), fill = "grey70") +
    geom_line(aes(x = inputs, y = mean)) +
    geom_point(data = df_refocus %>% dplyr::filter(trim2 == 0), aes(x = input, y = output), size = 5, shape = 1) +
    geom_point(data = df_refocus %>% dplyr::filter(trim2 == 1), aes(x = input, y = output), size = 5, shape = 4, color = "red4") +
    scale_y_continuous(breaks = seq(-1,5,1), limits = c(-2, 5)) + labs(y = "outputs") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 10))
# Final plot.
gridExtra::grid.arrange(p_refocus_1, p_refocus_2, p_refocus_3, ncol = 3)
```

<img src="HistoryMatching_explainer_files/figure-gfm/refocussing-1.png" height="80%" style="display: block; margin: auto;" />
As expected, our green region of Not Ruled Out Yet input-output tuples
shrinks as we repeatedly trim the unacceptable inputs. You can also see
that the model we re-fit to the trimmed data (black line) is changing
and regions of uncertainty (grey regions) are growing in the portions of
the input space that we trimmed. <br/><br/> <br/><br/>

### Matching only where we need to match.

It might be the case that all you’ve been given some historical
observations to use for matching, in which case you just have to hope
that they cover the range of inputs and outputs that you are interested
in. But given that you are trying to reduce uncertainty in your model’s
outputs, it’d be wise not to just randomly match input-output tuples for
history matching. It would be better to intentionally match input-output
tuples in areas of high output uncertainty in the hope of reducing the
uncertainty in the model output. (You won’t get rid of it all because
there is natural variance in underlying phenomenon. If it fits
perfectly, it’s very likely wrong.) <br/><br/> For example, I might fit
a Gaussian process emulator to the relationship between an input and an
output, which provides me with a cloud of uncertainty around an average
function that satisfies the observations:

``` r
# Generate some fake data.
inputs <- seq(from = 0, to = 20, by = 4)
inputs <- inputs[c(1:4,6)]
outputs <- sin(inputs * 0.4)
df_sim <- data.frame(inputs, outputs)
# Fit a Gaussian process emulator.
mod_GPE <-
  RobustGaSP::rgasp(design = inputs,
                    response = outputs)
predict_input <-
  data.frame( inputs = seq(from = min(inputs), to = max(inputs), by = 0.01))
df_GP <-
  data.frame(
      predict_input = predict_input,
      model_prediction = predict(mod_GPE, predict_input)
    )
# Plot.
(p_GPE <-
  ggplot() +
    geom_ribbon(data = df_GP, aes(x = inputs, ymin = model_prediction.lower95, ymax = model_prediction.upper95), fill = "grey70") +
    geom_point(data = df_sim, aes(x = inputs, y = outputs), size = 5, shape = 1) +
    geom_line(data = df_GP, aes(x = inputs, y = model_prediction.mean), linewidth = 1) +
    labs(title = "Matched points and emulated line\nwith 95% confidence intervals") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_blank(),
          text = element_text(size = 20))
)
```

<img src="HistoryMatching_explainer_files/figure-gfm/GPE-1.png" width="50%" style="display: block; margin: auto;" />
There are two things I want you take away from this plot:

1.  Note that the clouds of uncertainty disappear at the points where we
    have observed the input-output relationship, but expand in the
    regions between.

2.  Note that the cloud of uncertainty around the region of
    $12\le x \le 20$ is larger than elsewhere.

Let’s assume that we want to minimise uncertainty about the output
across the entire range of inputs. The plot above suggests that it would
be wise to observe outputs from inputs in the $12\le x \le 20$ range
because that is where our understanding of the outputs is worst. Once
observed, these would be *historical* outputs and ideal for history
matching. Alternatively, I can imagine a situation where our input
variable is time and we are more concerned about the near future than
the far future. In this scenario, we might choose to further minimise
the uncertainty about the output in the left portion of the plot and
worry less about the right portion.

``` r
# Update the fake data for the two scenarios.
inputs_early <- c(inputs, 1.5, 3.5, 6) %>% sort()
inputs_late <- c(inputs, 16) %>% sort()
outputs_early <- sin(inputs_early * 0.4)
outputs_late <- sin(inputs_late * 0.4)
df_sim_early <- data.frame(inputs = inputs_early, outputs = outputs_early, group = "Early matching")
df_sim_late <- data.frame(inputs = inputs_late, outputs = outputs_late, group = "Late matching")
df_sim_combined <-
  dplyr::bind_rows(df_sim_early, df_sim_late)
# Fit Gaussian process emulators to each.
mod_GPE_early <-
  RobustGaSP::rgasp(design = inputs_early,
                    response = outputs_early)
predict_input_early <-
  data.frame(inputs = seq(from = min(inputs_early), to = max(inputs_early), by = 0.01))
df_GPE_early <-
  data.frame(
      predict_input = predict_input_early,
      model_prediction = predict(mod_GPE_early, predict_input_early),
      group = "Early matching"
    )
mod_GPE_late <-
  RobustGaSP::rgasp(design = inputs_late,
                    response = outputs_late)
predict_input_late <-
  data.frame(inputs = seq(from = min(inputs_late), to = max(inputs_late), by = 0.01))
df_GPE_late <-
  data.frame(
      predict_input = predict_input_late,
      model_prediction = predict(mod_GPE_late, predict_input_late),
      group = "Late matching"
    )
# Plot.
df_GPE_combined <-
  dplyr::bind_rows(df_GPE_early, df_GPE_late)
(p_choosing_ctrl_pts <-
  df_GPE_combined %>%
  ggplot() +
    geom_ribbon(aes(x = inputs, ymin = model_prediction.lower95, ymax = model_prediction.upper95), fill = "grey70") +
    geom_point(data = df_sim_combined, aes(x = inputs, y = outputs), size = 5, shape = 1) +
    geom_line(aes(x = inputs, y = model_prediction.mean), linewidth = 1) +
    labs(x = "Input value", y = "Output value",
         title = "Choice of matching points and emulated line with 95% confidence intervals") +
    facet_wrap(~group) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 10))
)
```

<img src="HistoryMatching_explainer_files/figure-gfm/Two choices of adding matching points-1.png" style="display: block; margin: auto;" />
As you can see from the plots immediately above, sometimes we can choose
what to observe for our history matching in order to optimse particular
regions of a model’s performance. <br/><br/> <br/><br/>

### Continuous rather than batch updating

With methods like the implausibility scoring and the Not Ruled Out Yet
space, we have to run our model many times to map (cartographically
rather than mathematically) the relationship between inputs and outputs.
This can easily mean we do many unnecessary calculations. It is also a
retrospective method, which means we might spend a substantial amount of
time / iterations, not optimising our performance. We’d improved in
sudden jumps after periods of no improvement. Alternatively, we could
try to gradually optimise our performance something more akin to
real-time updating of our inputs, continuous rather than batch
processing. <br/><br/> So, what methods are available for continuously
updating our inputs? <br/><br/>

#### Bayesian calibration

Bayesian calibration involves the application of Bayes’ rule, which is
explained in a million places online and those good ol’ fashioned things
called books. Very briefly, we start with some initial guess at the
distribution of input values that result in acceptable outputs, and then
update that distribution with observations of inputs that we know to
result in acceptable outputs. <br/><br/> One distinction between typical
history matching and Bayesian calibration is somewhat analogous to the
distinction between fixed effects and random effects in regression
analysis. The analogy I’m alluding to is that estimating regression
coefficients for fixed effects is like saying “*This is the estimate;
other estimates are not acceptable*”, and filtering parameter values
using history matching is like saying “*These parameter values are
acceptable. Other parameter values are not acceptable*”. <br/><br/>
Similarly, estimating regression coefficients for random effects is like
saying “*This is the distribution of acceptable coefficients, but I’m
not saying any particular one is best*”, and filtering parameter values
using Bayesian calibration is like saying “*This is the distribution of
acceptable parameter values, but I’m not saying any particular one is
best*”. <br/><br/> (*Side note: Of course, there are many differences
between fitting regression coefficients and filtering parameter values
in simulation models. I present the previous analogy to share something
that helped me on my journey to understand all these concepts. Analogies
are temporary rafts to help us across rivers of confusion; they need not
be perfect and they can be left behind once we reach the firm
understanding of the other side.*) <br/><br/> Bayesian calibration, can,
in principle, handle what history matching can’t by dealing with
distributions rather than examples. But it requires us to have done a
good job of incorporating the likelihood of the aforementioned outputs
in the prior. This distinction between history matching and Bayesian
calibration might be pithily summarised as “History matching is an
easier way to learn about a little, while Bayesian calibration is a
more-difficult way to learn about a lot”. <br/><br/> [Edwards, Cameron
and Rougier](https://www.sci-hub.wf/10.1007/s00382-010-0921-0) suggest
that history matching can be a useful pre-calibration step to sensibly
inform the prior for Bayesian calibration, but this is only reasonable
if the desirable responses of the outputs to the inputs is centred
within the Not Yet Ruled Out space. By analogy, sure, it seems sensible
to start looking for your lost keys under the streetlight at night, but
that doesn’t mean they aren’t *actually* somewhere in the dark.
<br/><br/> [Williamson et
al.](https://www.sci-hub.wf/10.1007/s00382-013-1896-4) provides a small
discussion on differences between history matching and Bayesian
calibration, if you are interested. <br/><br/> <br/><br/>

#### Ensemble Kalman filters

But perhaps you are more of a frequentist than a Bayesian. Fear not;
there is a gradual-update approach for you, too. Well, actually, it
still uses Bayes theorem but you don’t have to specify the form of the
prior, which is the difficult part. Instead, you just assume it is
Gaussian and that it is defined by the average and variance of some
previous, small searches. <br/><br/> Ensemble Kalman filtering involves
iterating two-steps to update our behaviour in response to our
observations. THis kind of filtering is suited to systems where the
outcome influences the subsequent input cyclically, rather than the
input merely nudging an on-going output-generating process. For example,
have a read of [Ward, Evans, and Malleson
2016](https://sci-hub.wf/10.1098/rsos.150703) to see how ensemble Kalman
filters are used to dynamically calibrate agent-based models. <br/><br/>
There are two steps:

1.  forecasting

2.  assimilating observations

The first step is to calculate the immediately-future state of your
model, for a variety of possible inputs. It’s like when you just heard
your phone buzz from receiving a message but you can’t recall where you
left it. At first, you sit still and limit yourself to observing only
few feet around you in many directions. Did the buzz come from the left
a little, to the right a little, etc.? This gives a set of possible near
next moves - a.k.a. an ensemble. Proponents of ensemble Kalman filters
suggest that the mean of this ensemble is the best guess of the desired
output, and the variance is a measure of uncertainty. <br/><br/> The
second step is to actually observe the immediate future. This is like
your phone buzzing again as a second message arrives, but you are paying
attention this time so you now know for sure that it is somewhere to the
right of you. Congratulations, you have updated your understanding and
limited the range of possible inputs toward your desired output. In
other words, you have constrained the choice of directions in which you
will travel to find your phone. <br/><br/> The fossil fuel industry seem
to be rather keen on this approach - see [Rwechungura et
al. 2011](https://sci-hub.wf/10.2118/142497-MS). I’m unsettled about the
Gaussian assumption but, as Anna said in Frozen 2…

``` r
knitr::include_graphics("https://media.giphy.com/media/Ymt6N7O93ixVVbmBNl/giphy.gif")
```

<img src="https://media.giphy.com/media/Ymt6N7O93ixVVbmBNl/giphy.gif" width="50%" style="display: block; margin: auto;" />

I have to agree with her: in times when you can’t know the future with
sufficient certainty, choose the least-worst option, take a step, then
reassess. <br/><br/> (*Side note: This is exactly the approach advocated
for navigating complex adaptive systems, particular social ones. [Dave
Snowden](https://thecynefin.co/team/dave-snowden/), coiner of the term
“antropic complexity”, [is particularly a
fan](https://www.coachesrising.com/podcast/the-problem-with-developmental-models-with-dave-snowden/).*)
<br/><br/> <br/><br/>

## Final thoughts and other resources

Hopefully, you now know what history matching is and why it is useful.
One irksome thought I had while reading up on history matching was
*Won’t it be really expensive to keep running loads of iterations of the
model in order to optimise the inputs?* If you thought the same, then
rest assured that we are not alone. This concern has led people to use
emulators to reduce the computational burden. These are models of
models, which you can read more about in my other blog about [Gaussian
process emulators](https://github.com/ciaranmci/SoMaS). <br/><br/> Below
are my thoughts that I’d like to discuss with a history-matching expert,
someday.

### The dangers of selection

History matching assumes that previous outputs are sufficient proxy
observations of the truth, so we should use them to constrain our wider
understanding of the truth. But what if our previous observations of
outputs poorly represent the truth? Selection bias is a beast!
<br/><br/> I think selection bias is easiest to understand as collider
bias using directed acyclic graphs, which is beyond the scope of this
blog. In place of a deeper dive, check out figure 1 in [Griffith et
al.](https://www.sci-hub.wf/10.1038/s41467-020-19478-2) for a simple
illustration of how selection bias changes our perception of a
relationship (image URL is in code box). Focus on part C of the figure
below:

``` r
knitr::include_graphics("https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-020-19478-2/MediaObjects/41467_2020_19478_Fig1_HTML.png?as=webp")
```

<img src="https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-020-19478-2/MediaObjects/41467_2020_19478_Fig1_HTML.png?as=webp" alt="Image is figure 1 taken from [Griffith et al. (2020)](https://sci-hub.wf/10.1038/s41467-020-19478-2) '_Collider bias undermines our understanding of COVID-19 disease risk and severity_'" style="display: block; margin: auto;" />
<br/><br/>

[Williamson et al.](https://www.sci-hub.wf/10.1007/s00382-013-1896-4)
hint at something like selection bias when they say the following about
ensemble filtering:

> “The issue with this approach is that it does not account for the
> uncertainty in the unsampled, high dimensional parameter space.”

This uncertainty means that our historical observations are random
samples of an uncertain variable, yet we use the historical observation
as if it were not. This is akin to the way we assume our covariates /
independent variables in regression modelling are measured without
error, which is obviously unjustifiable (*Side note: Regression models
that permit error in the independent variables [do
exist](https://en.wikipedia.org/wiki/Errors-in-variables_models) but
unfortunately are rarely used.*) <br/><br/> To be fair, the issue of
sampling from an uncertain variable is brushed away simply by repeated
sampling and taking averages. This is different to selection bias, which
typically refers to bias arising from focusing on only a stratum of a
variable’s values. Averaging won’t alleviate the resulting bias because
the sampling is not completely random. But if we only take a few
observations for our history matching (rather than an average of
repeated samples), then we are committing the crime of selection bias by
effectively only sampling from a portion of the variable’s values. I’m
glad Williamson et al. say what they do about ensemble filtering but
it’s a shame they don’t acknowledge this issue for history matching; why
not? Perhaps I am missing something. Please, let me know. <br/><br/>
*(Side note: Focussing our analyses is permissible if we acknowledge
that we are doing it. But if you don’t realise that you have selected a
subset of situations and you make statements about the whole, then your
inferences are incorrect. What I’m trying to say is that bias is not
inherently a bad thing, but not acknowledging it misleads.)* <br/><br/>

### The dangers of extrapolation

The models we fit usign history-matched points make “predictions” beyond
the matched points, i.e. interpolation and extrapolation. Another
problem I have with history matching is that it infers what inputs are
reasonable by assessing only a handful of observations. Making
inferences is just part of life but it must be handled with care. For
example, does history matching assume a linearity between input and
output values? By disregarding all input values beyond our Not Ruled Out
Yet space, we don’t permit non-linear effects, reversal of effects, or
pockets of sensible outputs after a gap of nonsensical ones. Let me try
to give some illustrative examples:

-   An over-the-top example is to match a straight-line linear
    regression model to pre-1960 global surface temperature. This is
    something like history matching on only these pre-1960 data. This
    underestimates today’s temperatures. This could be avoided by not
    limiting our matching observations to a small region. Latin
    Hypercube sampling and theoretically-informed ranges should be able
    to account for this, in many cases. But what if all historical
    observations - and thus our understanding - manifest from the same
    local portion of a curve?

``` r
# Prepare data
webscrape <- read.csv("https://data.giss.nasa.gov/gistemp/graphs/graph_data/Global_Mean_Estimates_based_on_Land_and_Ocean_Data/graph.txt")
df <-
  webscrape %>% 
  tidyr::separate_rows(colnames(webscrape), sep = ", ") %>% 
  dplyr::slice(4:n()) %>%
  tidyr::separate(colnames(webscrape), sep = "\\s+", into = c("Year", "No_Smoothing", "Lowess(5)")) %>%
  dplyr::select(-`Lowess(5)`) %>%
  dplyr::mutate_if(is.character, as.numeric)
# Add the 1951-1980 baseline value
# https://earthobservatory.nasa.gov/world-of-change/decadaltemp.php
df$No_Smoothing <- df$No_Smoothing + 14
# Fit linear linear model to data from 1960.
timespan <-
  df %>% dplyr::select(Year)
df <-
  df %>%
  dplyr::filter(Year < 1920) %>%
  lm(formula = No_Smoothing ~ Year, data = .) %>%
  predict(newdata = df %>% dplyr::select(Year)) %>%
  as.data.frame() %>%
  `colnames<-`("trend") %>%
  dplyr::bind_cols(df)
# Plot.
(p_tempExtrap <-
    df %>%
    ggplot() +
     geom_point(aes(x = Year, y = No_Smoothing),
                size = 3, color = "grey") +
     geom_line(aes(x = Year, y = trend),
               linewidth = 3) +
    ylab("Global average\nsurface temperature (Celsius)") +
    labs(caption = "Data source: https://climate.nasa.gov/vital-signs/global-temperature/") +
    theme(text = element_text(size = 10))
    
)
```

<img src="HistoryMatching_explainer_files/figure-gfm/Global temp extrapolation -1.png" width="50%" style="display: block; margin: auto;" />

-   Conor Crilly beautifully shows wide uncertainty at both ends of an
    emulated signal, in his [2022 blog about Gaussian process
    emulation](https://compass.blogs.bristol.ac.uk/2022/01/25/gaussian-process-emulation/).
    This could be avoided by constraining the emulation, as suggested by
    Crilly.

``` r
knitr::include_graphics("https://bpb-eu-w2.wpmucdn.com/blogs.bristol.ac.uk/dist/e/692/files/2022/01/PlotTwo-1024x512.png")
```

<img src="https://bpb-eu-w2.wpmucdn.com/blogs.bristol.ac.uk/dist/e/692/files/2022/01/PlotTwo-1024x512.png" style="display: block; margin: auto;" />

-   [Johnny
    Hyman](https://github.com/jonnyhyman/Chaos/blob/master/logistic_zoom.py)’s
    simulation of logistic map shows how unpredictability arises from a
    simple rule: $X_{n+1} = r \times X (1 - X)$. In other words,
    [deterministic chaos](https://en.wikipedia.org/wiki/Chaos_theory).
    There is no way that historical observations can help us, here. As
    many people have said, it’s like driving forward while looking in
    the rearview mirror. This simulation was made for a [Versatium
    video](https://www.youtube.com/watch?v=ovJcsL7vyrk&t=766s) that is
    well worth watching.

``` r
knitr::include_graphics("https://i.makeagif.com/media/6-28-2023/O09D6n.gif")
```

<img src="https://i.makeagif.com/media/6-28-2023/O09D6n.gif" width="60%" style="display: block; margin: auto;" />
<br/><br/>

### Baby out with the bathwater

Finally, as a parting idea…the gist of [Williamson et
al.](https://www.sci-hub.wf/10.1007/s00382-013-1896-4)’s idea of the Not
Ruled Out Yet space is that we throw out the input values that *can*
produce ridiculous and unlikely outputs. But some of these input values
might have produced output values within the Not Ruled Out Yet space. It
seems very conservative to throw out any input values that misbehaved a
little. <br/><br/> <br/><br/>

### Resources

1.  The function I used in the example of history matching with absolute
    values was taken from the tutorial for the [`hmer`
    package](https://danny-sc.github.io/Tutorial_1/one-dimensional-example.html),
    which is a history-matching package in R. If you want a Python
    package for history matching, check out
    [`mogp_emulator`](https://mogp-emulator.readthedocs.io/en/latest/demos/historymatch_demos.html).
