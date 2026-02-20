# CellRank Protocol Teaching Blueprint

This blueprint turns the project into a complete teaching track for:

`03_CellRank_protool.pdf`  
(`/Users/ashi/Library/CloudStorage/OneDrive-UAB-TheUniversityofAlabamaatBirmingham/0Informatics Club/CCC/03_CellRank_protool.pdf`)

It is designed to ensure all protocol sections are taught, not only code execution.

## Primary teaching notebook

Use the workshop notebook as the practical backbone:

- `../sctalk3_cellrank2/cellrank2_tutorial.ipynb`

## Coverage map (Protocol -> Lesson in this project)

| Protocol section (PDF) | Where it is taught | Required learner output |
|---|---|---|
| Introduction | Notebook: `Part 0` + opening sections | Explain why fate mapping needs Markov chains in single-cell data |
| Overview of procedure | Notebook: `Part 0` + `Summary & Key Takeaways` | Draw the 3-stage CellRank workflow (kernels -> estimator -> downstream) |
| Input requirements | Notebook: each procedure “About the Dataset” + setup | Identify required `AnnData` slots (`obs`, `obsm`, `obsp`, `layers`) per kernel |
| Stage 1: kernels | Procedures 1-4 | Build transition matrices with `CytoTRACEKernel`, `PseudotimeKernel`, `VelocityKernel`, `RealTimeKernel` |
| Optional kernel combination | Procedure 3 | Combine kernels (e.g., velocity + connectivity), justify weighting choice |
| Stage 2: estimator (GPCCA) | Procedures 1-4, steps for macrostates/terminal states | Compute macrostates, initial/terminal states, and defend chosen `n_states` |
| Stage 3: downstream analyses | Procedures 2-4 | Compute fate probabilities, lineage drivers, and gene trends |
| Limitations | This blueprint (section below) + in-class discussion | Submit a limitation/risk memo for one dataset |
| Experimental design | This blueprint (section below) | Draft a kernel-selection rationale before running code |
| Materials/Equipment/Software | `environment.yml`, setup sections, conda CLI | Reproduce environment and run all imports without error |
| Troubleshooting | This blueprint (section below) | Resolve at least one induced failure mode and document the fix |
| Timing | This blueprint (section below) | Produce a realistic run plan for a workshop session |
| Anticipated results | This blueprint + notebook outputs | Deliver expected plots/tables and interpret biological meaning |
| Data availability / Code availability | Notebook references + repo files | Locate and cite dataset + code sources for each procedure |

## Instructor flow (teach all protocol content)

1. Concept foundation
   - Teach Markov chain intuition, transition matrices, macrostates, absorbing states, and fate probabilities.
2. Stage 1 in practice
   - Run each kernel on its matched dataset and compare assumptions.
3. Stage 2 in practice
   - Use GPCCA to compute macrostates, terminal states, and initial states.
4. Stage 3 in practice
   - Compute lineage drivers and gene trends, then interpret biological trajectories.
5. Protocol-critical non-code topics
   - Limitations, troubleshooting strategy, timing, and experimental design trade-offs.
6. Validation and rigor
   - For LARRY, compare inferred fate with lineage tracing outcomes.

## Limitations to teach explicitly

- RNA velocity assumptions can fail on low-quality splicing signal or weak dynamical fit.
- Pseudotime can impose a misleading direction when branching complexity is high.
- OT coupling quality depends on time-point coverage and preprocessing choices.
- Terminal state calls depend on macrostate resolution (`n_states`) and may be unstable.
- Driver ranking is associative, not causal; follow-up validation is required.

## Experimental design checklist (before running)

- Define biological question and expected terminal lineages.
- Select kernel(s) based on available data views:
  - velocity present -> `VelocityKernel`
  - pseudotime available -> `PseudotimeKernel`
  - time-resolved data -> `RealTimeKernel`
  - no directional prior -> `CytoTRACEKernel`
- Define validation strategy:
  - known markers, known transitions, lineage tracing, or external benchmarks.
- Predefine sensitivity checks:
  - vary `n_states`, neighbors graph parameters, and kernel combinations.

## Troubleshooting drills (must be taught)

- Missing graph in `obsp` -> rebuild neighbors/connectivity graph.
- Missing pseudotime key in `obs` -> recompute or rename key consistently.
- GPCCA decomposition instability -> change `n_states`, inspect coarse-graining quality.
- Empty or weak driver output -> verify terminal states, lineage labels, and expression preprocessing.
- OT/RealTime issues -> verify time labels ordering and moscot coupling construction.

## Timing template (single workshop day)

- Module 1 (30 min): Concepts + architecture
- Module 2 (45 min): Procedure 1 (CytoTRACE)
- Module 3 (45 min): Procedure 2 (Pseudotime + drivers + trends)
- Module 4 (45 min): Procedure 3 (Velocity + kernel combination)
- Module 5 (45 min): Procedure 4 (RealTime + OT + validation)
- Module 6 (30 min): Limitations, troubleshooting, design review

## Completion rubric (all required for “teach all of these”)

- [ ] All four procedures executed end-to-end at least once.
- [ ] Learner can justify kernel choice for a new dataset.
- [ ] Learner can interpret macrostates, terminal states, and fate probabilities.
- [ ] Learner can produce and interpret lineage driver genes and gene trends.
- [ ] Learner can describe at least three protocol limitations and mitigations.
- [ ] Learner can troubleshoot at least one failure mode independently.
- [ ] Learner can report reproducible runtime plan and environment details.

## Minimal run commands

```bash
conda env update -f environment.yml --prune
conda activate cellrank2exp
cellrank2exp start
```

For notebook execution, run Jupyter from the repo root and open:

- `sctalk3_cellrank2/cellrank2_tutorial.ipynb`
