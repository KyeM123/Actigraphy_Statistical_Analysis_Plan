# Actigraphy_Statistical_Analysis_Plan
The purpose of this script is to address sleep data collected via wristwatch actigraphy. 
Authors: K.M. and R.V.R.
Date: 28/05/2026

The purpose of this script is to address sleep data collected via wristwatch actigraphy. 

The study is a double-blinded placebo-controlled crossover RCT, in which participants were allocated to one of two groups. Each group received a placebo or active treatment for six weeks, before undergoing a four-week washout period and receiving the other treatment for a further six weeks. Wristwatch actigraphy measures were collected for the seven days preceding each trial arm and for the last seven days of each arm. Survey data were also collected pre and post trial arms.

The primary outcomes were objective actigraphy-derived sleep parameters, specifically sleep period time (SPT), sleep duration, sleep efficiency and wake after sleep onset (WASO). These measures capture complementary aspects of sleep quantity and quality and are therefore analysed concurrently. Secondary outcomes include subjective sleep quality as assessed using the Pittsburgh Sleep Quality Index (PSQI). Secondary outcomes will be analysed in a separate script.

**Statistical Analysis**
The primary objective of the analysis is to assess whether pre-to-post changes differ between treatment conditions, expressed as a difference-in-differences (DID) effect. 

The dataset contains repeated observations for each participant across:
-	Treatment (within-subject): A vs B
-	Time (within-subject): pre vs post
-	Treatment Period (within-subject): period 1 vs period 2
-	Treatment Order (between subject): A-first vs B-first

For each outcome, linear mixed effects models with a random intercept (1|id) are fitted. Fixed effects structures, including treatment x time with optional adjustments for order and period are compared using AICc, with the best-fitting model for each outcome selected. The models will then be tested for non-normality, heteroscedasticity and collinearity. If the model assumptions are violated, transformations are applied sequentially (square-root then log), followed by a gamma generalised linear mixed model. Assumptions are also evaluated manually through visual inspection, as formal statistical tests can be overly sensitive with small sample sizes and sleep data. There is an option to override the automatic model selection based on visual diagnostics.

The primary treatment effect is assessed via the treatment x time interaction from Type III ANOVA. DID estimates and 95% confidence intervals are obtained using estimated marginal means.

Sensitivity analyses are conducted using Cook’s distance based on leave-one-ID-out (LOIO) refitting. For flagged participants, LOIO analyses are performed to assess the robustness of DID estimates.
