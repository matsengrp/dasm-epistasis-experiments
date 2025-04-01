# Project overview
Use DASM to identify and calculate epistasis between positions and mutations in antibodies.

## What are you trying to do? Articulate your objectives using words that would be familiar to someone who had taken an undergraduate class in the topic.
We aim to identify epistasis between positions/mutations within a single antibody sequence using predictions from the DASM model (which predicts selection factors for every position in an antibody according to sequence). Epistasis occurs when the effect of a mutation at one position depends on the amino acids present at other positions. Our objective is to systematically quantify these non-additive interactions in the antibody sequence to understand how different positions functionally interact with each other. By calculating epistasis values between positions and analyzing their statistical significance, we aim to uncover networks of functionally coupled residues that contribute to antibody stability, binding affinity, and specificity. This will provide insights into the structural and functional constraints on antibody evolution during affinity maturation.

## Who cares? If you're successful, what difference will it make?
This work can:
1. Reveal mechanisms driving antibody affinity maturation and evolution. Epistasis creates "constrained pathways" where certain mutations only become beneficial after others have occurred, creating a sequential pattern of mutation accumulation that couldn't happen through purely additive effects.
2. Reveal structural and functional insights and constrictions about antibodies, assisting in clinical development of antibody engineering and vaccine development.

## How is it done today, and what are the limits of current practice?
Current approaches to studying antibody epistasis include:
1. Deep mutational scanning (DMS) experiments that measure binding or expression for thousands of variants. Limitation: are labor-intensive and so limited to specific antibodies and antigens.
%EM Also, most DMSs are point mutations only, right?
2. Computational predictions based on structural models are possible, but these might not be accurate enough to catch the complex interactions within the protein in its bound and unbound states.
%EM In any case it seems like evolutionary information is a complement to structural information. Hugh pointed out the other day that things that work well to stabilize a structure (e.g. hydrophobics on the outside) can cause bad biophysical properties in the real world.

## What's new in your approach and why do you think it will be successful?
It uses a model (DASM) that is trained on antibody sequences and today does the most accurate job at predicting selection factors per sequence.
We rely on this model to test the relationship between pairs of mutations and the effects of a mutation on different backgrounds.

## What is a best case scenario hypothetical result? Try to be as specific as possible. Use your imagination!
We are able to predict epistasis in any antibody sequence, and our predictions are experimentally validated as correct.
%EM We know Johanne Jacobsen https://www.ous-research.no/jjacobsenlab who seems willing to make and test antibodies for us. Also David Glass (see email forward).
Vaccine and theraputic antibody design are shifted from a "one mutation at a time" approach to a more sophisticated view that leverages the complex dependencies within antibody structure and function, informing strategies such as sequential therapeutic interventions

## What are the potential bad outcomes? Any overall concerns here?
- We are more influenced by noise in the data/not enough data than by the actual biology.
- Uncertainty measurements are too high to rely on the predictions.
- Our predictions do not match existing experimental inferences.

## Is there pre-existing work/code that could be leveraged to explore the potential for bad outcomes? To do proof-of-concept investigation to get a first-pass answer for the underlying scientific question?
The project is based upon the DASM model predictions. For a proof-of-concept, once simple predictions can be made, we should start by testing antibody sequences with known epistatic relationships, for example:
- https://www.pnas.org/doi/10.1073/pnas.2413884122, Schulz 2025 - experimentally test 2^10=1024 variants (variable in the heavy chain only) of a SARS-CoV-2-specific antibody, COV107-23, that targets the receptor binding domain of the viral spike protein, and find an epistasis hotspot at residue 53 (also an interesting paper as they fit a pairwise and global model to their experimental results).
- https://elifesciences.org/articles/83628.pdf, Phillips 2023 - experimentally test antibody library of CH65 antibody, which is a broad antibody to diverse H1 influenza strain. They identify epistasis both in the heavy chain, light chain, and between both. This group has other papers about antibody epistasis as well.

%EM These are great. It would be neat if we could train on some covid repertoires for the first (we have some, called Ye) and some flu repertoires for the second (ask Mackenzie?).

## Are there any other categorically different approaches that could be applied here?
Not that I am aware of.
%EM Well, if one had a structure one could use https://www.science.org/doi/10.1126/science.add2187 or https://www.biorxiv.org/content/10.1101/2024.07.09.602403v2 . It would be interesting to compare to these, both with a "real" PDB structure and also with an AlphaFold structure. The authors of the Hermes, the latter, are "lab friends" and are keen to collaborate. Hugh knows the MPNN folks.

%EM More about Hermes: https://www.pnas.org/doi/10.1073/pnas.2300838121 . The best intro is probably her video https://online.kitp.ucsb.edu/online/viralimmune24/nourmohammad/rm/jwvideo.html 

## If this is a methods project, what methods will you compare to? Can you get them running before writing new code?
Not a methods paper (?)
%EM Above

## What data will you use? Are there appropriate hold-out sets?
Predictions will be made using existing DASM models. 

Validations for pairwise epistasis:
- experimental data such as mentioned three questions back
- structurally resolved antibodies
- co-occurence and *order of occurence* of mutations in existing antibody trees.

%EM I wonder how Hermes would do on the Schulz data. We could start a chat with them.

Validations for background dependent epistasis:
- data from Whitehead lab that does DMS on 9 different genetic backgrounds (https://www.nature.com/articles/s41467-024-48072-z)
- DMS data on different antibodies 

## Is it possible that better data would make this project irrelevant?
Not soon:
- Pairwise epistasis experimental measurements for a such a large antibody space are not possible with the current technologies. 
- Specific interactions and cases can be tested experimentally as a more accurate test, but these specific probes can be directed by a computational method such as the one proposed here.

## Sketch the approach, broken down into steps, with expected amounts of time and intermediate steps for each.

### Stage 1 (1 months):
Write code that calculates epistasis by perturbing positions. Specifically:
- Take a sequence. 
- For every position in the sequence, mutate it into every possible amino acid and measure deviations in the other positions. 
- Calculate some measure of uncertainty for these predictions. 
- Positions that change the predictions in other positions in unexpected ways could be epistatic hot spots. 

Use method on known sequences with experimentally proven epistatic relationships. Does this work? (POC)
- If yes, move to next step.
- If not, stop and rethink.

### Stage 2 (1 month)
- Create similar analysis for background dependent epistasis.
- Further validations on experimental data for both on datasets.
  - Obtain and/or organize Whitehead lab data.

### Stage 3 (?)
Calculate co-occurence and order of occurence for mutations in antibody trees. 
See if epistatic patterns found in DASM match this.

### Stage 3 (?)
If we have systemic predictions of epistasis that work, we can use this to widen our understanding of their general evolutionary process:
- How prevalent is epistasis and in what sections of the antibody?
- What order interactions do we see (first, >2)
- Are the epistasis interaction global or local? See https://www.pnas.org/doi/full/10.1073/pnas.1804015115 and https://www.pnas.org/doi/10.1073/pnas.2413884122. 

#### How do we decide to move onto the next stage?

