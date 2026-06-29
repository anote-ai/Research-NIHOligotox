# Can AI predict whether a gene-silencing drug will be toxic -- before it's ever made?

*A plain-language look at the OligoTox project*

## The problem: a promising class of drugs keeps failing on safety, not effectiveness

Over the last decade, a new class of medicines called **oligonucleotide
therapeutics** -- short synthetic snippets of DNA or RNA -- has moved from lab
curiosity to approved drug. They work by directly silencing or correcting the
genetic instructions behind a disease, rather than blocking a protein the way most
traditional pills do. More than 15 of these drugs are approved, and hundreds more
are in clinical trials, for everything from spinal muscular atrophy to high
cholesterol.

The problem is that many of these drug candidates don't fail because they don't
work -- they fail because they're **toxic**. The liver and kidneys can react badly
to the chemical backbone of these molecules. The immune system can mistake certain
DNA sequences for an infection and trigger a dangerous inflammatory response.
Platelets can clump or the blood can stop clotting properly. Today, most of this is
discovered the expensive, slow way: in animal studies, late in the drug development
pipeline, after a lot of money has already been spent.

## The idea: teach a model to "read" both the chemistry and the sequence

Existing computer models that try to predict this kind of toxicity mostly look at
the drug's chemistry -- its molecular fingerprint -- and ignore something equally
important: the actual genetic sequence of letters (A, T, G, C) that makes up the
drug. That sequence matters because it determines what *else* in the body the drug
might accidentally stick to (called "off-target binding"), and because certain
sequence patterns (like repeated "CG" pairs) are known to set off immune alarms.

This project's core scientific bet, laid out in `DESIGN_DOC.md`, is that combining
**sequence** information with **chemistry** information -- an approach called
"SeqChem" -- should predict toxicity meaningfully better than chemistry alone,
especially for genuinely novel drug designs the model hasn't seen anything similar
to before (the hardest, most realistic test of whether a model actually generalizes,
versus just memorizing).

## Where the project actually stands today

It's worth being candid about what exists right now versus what's still on the
drawing board, because this is an active research project, not a finished result.

**What's built and running:**
- A small, working Python toolkit (`src/oligotox/`) that can take a DNA/RNA
  sequence, compute basic features (GC content, CpG content, length, backbone
  chemistry), and feed them into a simple classifier.
- A synthetic data generator that creates fake-but-structured example
  oligonucleotides for testing the pipeline end-to-end.
- A baseline experiment (`experiments/exp0_baseline.py`) that actually runs this
  pipeline and reports a real, measured accuracy number on the synthetic data.
- A separate, larger effort: a proposal submitted to the **NIH NCATS OligoTox Open
  Data Challenge** (see `OligoTox_Phase2_Submission_Anote.md`) to generate a much
  bigger, *real*, human-cell-based toxicity dataset -- 2,800 drug-like molecules
  tested across 47 different toxicity measurements in five different lab systems.
  This project was a Phase 1 winner of that NIH challenge.

**What's still aspirational:**
- The actual "SeqChem" combined representation described above hasn't been built
  yet -- today's toolkit uses a handful of simple features, not the richer
  sequence-embedding-plus-chemistry-fingerprint design described in the research
  plan.
- No model has yet been trained or tested on real-world toxicity data. All current
  results come from synthetic, programmatically generated examples, which is useful
  for confirming the code works, but doesn't tell us anything about real biology yet.
- The "hardest test" mentioned above -- checking whether a model still works well on
  genuinely novel drug designs it wasn't trained on -- hasn't been implemented.

## Why this still matters, even unfinished

Even at this early stage, the project is doing something useful: building toward a
*public, shared benchmark* for oligonucleotide safety prediction, in a field where
most of the relevant safety data today is locked up inside individual pharma
companies. If the planned NIH-backed wet-lab dataset materializes on schedule (the
project's plan targets rolling data releases between August and November 2026), it
would become one of the largest open datasets of its kind, and would give the
broader research community something to train and test sequence-aware toxicity
models on for the first time.

The honest summary: the scientific hypothesis is well-motivated and grounded in
known biology (certain sequence motifs really do trigger immune responses; backbone
chemistry really does drive liver toxicity), the data-generation plan is credible
and already recognized by NIH, but the AI/ML side of the project -- the part that
would actually prove "sequence + chemistry beats chemistry alone" -- is still mostly
a plan rather than a result. The next milestone to watch for is the first real
experiment run on either the NIH wet-lab data or an interim public dataset, instead
of synthetic placeholders.

---
*This post summarizes the state of the `research-niholigotox` project as of this
writing. For technical details, see `DESIGN_DOC.md` (the AI/ML research plan) and
`OligoTox_Phase2_Submission_Anote.md` (the NIH data-generation submission).*
