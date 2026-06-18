# Agent Instructions for this project
 
## Project Overview
- This is an astronomy research project that analyzes the faithfulness of the Jackknife sampling method for galaxy cluster shear profile measurements, for a galaxy cluster cosmology analysis. The analysis results should be useful for Vera Rubin's LSST data. 

## Goals
- Goal 1: Derive mesurements of galaxy cluster shear profiles on simlations and their data-based covariance matrix measurements derived through jackknife measurements.
- Goal 2: Compare the variations of the measurements derived between simulations versus the jackknife-sampled covariance

## Technical Context
- **Language/Frameworks:** Python
- **Key Libraries:** treecorr(for profile measurement and covariance), astropy, matplotlib, numpy

## Coding Standards & Conventions
- **Style Guide:** Use functional components, be modular, and configurable
- **Best Practices:** Write test for specific functional component

## Constraints & Rules (Dos and Don'ts)
- **Do:** Write python scripts, and then execute the python scripts when doing analysis
- **Do:** Write your progress into a claude_progress.md file. 
- **Do:** Before ending each session, ask if you should commit all changes, and do so if given permission.
- **Do:** Do not git pull or push. Prompt the human developer to do so, if needed.
- **Don't:** Do not go outside the folder that Claude was started in.
- **Don't:** Do not submit batch jobs without asking.

## File Structure
- `configs/`: configuration files for specifiying the measurements
- `validations/`: validation analyses for testing the code
- `docs/`: Documentation

## Behavior Guidelines below:

### 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

### 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

### 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

### 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

---

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.
