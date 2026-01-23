
# scarscape development philosophy

`scarscape` is a **self-contained**, **streaming**, and **schema-stable** CLI for computing summary statistics from cancer-genome inputs (reference genome, small-variant VCF, SV VCF, copy-number segments, and region sets).

It is designed to be **boring to run**, **hard to misuse**, and **easy to extend**.

---

## Non-negotiables

### Self-contained binary
- The `scarscape` binary must run without system dependencies (no `bcftools`, no `tabix`, no shelling out).
- We may *recommend* upstream normalization (e.g. `bcftools norm`) for convenience, but `scarscape` never executes external programs.

### Streaming by default

Genome statistics do not require large amounts of memory to compute.

- `scarscape` should use negligable memory to leave room for parallelisation.
- Inputs are processed as streams/iterators wherever possible.
- Avoid loading whole VCFs or entire genomes into memory.

### Fail fast where correctness is threatened
- Invalid or inconsistent inputs should yield **actionable, sample-scoped errors** early.
- When an error compromises the correctness of a downstream statistic, the tool should stop (default).
- If an error only affects an optional metric, the tool may skip that metric and continue **only if** it can report the skip explicitly.

### Stable outputs and reproducibility
- Output table schemas (column names, meanings, types) are versioned and treated as a contract.
- Runs should be reproducible given the same inputs and reference bundle.

---

## The four concerns and their handoffs

`scarscape` is structured around four concerns. Each concern has a narrow responsibility and produces a well-defined artifact for the next layer.

### 1) I/O adapters (reading + validation + normalization)
**Responsibility**
- Read input files (manifest, VCFs, segments, FASTA, region sets).
- Perform format validation (e.g. required fields, monotonic constraints, basic invariants).
- Enforce “strictness” rules (e.g. reference compatibility checks).
- Apply **in-process normalization** where feasible and deterministic.

**Must not**
- Compute summary statistics.
- Know about output schemas.
- Accumulate large intermediate datasets.

**Handoff artifact**
- **Normalized event streams** (iterators) that yield validated domain events.
- Reference handles (indexed FASTA readers, contig maps, region set handles).
- Clear, structured errors that include: sample ID, file path, record location (when available).

> Contract: If the I/O layer yields an event, it is structurally valid and normalized to the degree required by downstream computations.

---

### 2) Domain model (meaningful types + invariants)
**Responsibility**
- Define the conceptual vocabulary of the program: sample identity, genome build, coordinates, event categories, region sets, and “what a normalized variant means”.
- Encode invariants in types where possible.
- Provide small, deterministic transformations that do not require external I/O.

**Must not**
- Read files directly.
- Write output tables.
- Contain “business logic” for any particular analysis.

**Handoff artifact**
- A set of domain types used across the codebase so higher layers do not pass around raw strings/maps.

> Contract: All layers agree on a single definition of coordinates, contigs, and normalized event semantics.

---

### 3) Stats engines (pure computations)
**Responsibility**
- Consume normalized event streams + reference metadata and compute summary statistics.
- Be deterministic, testable, and side-effect free.
- Prefer composable primitives (e.g. a reusable “count-in-region-set” engine used by many statistics).

**Must not**
- Re-parse raw file formats (no direct VCF field parsing here).
- Contain output formatting or file writing.
- Perform external I/O beyond what is required for reference lookups (e.g. querying indexed FASTA).

**Handoff artifact**
- Typed summaries (per-sample or per-run) suitable for serialization.
- QC metadata about what was computed, skipped, or filtered (when applicable).

> Contract: Stats engines assume inputs are already validated/normalized. This assumption should be guaranteed by the types leveraged in the adapter layer.

---

### 4) Reporting (schemas + serialization + file layout)
**Responsibility**
- Define stable table schemas and write outputs consistently.
- Handle formatting (CSV/TSV), ordering, headers, and version tags.
- Own the “what files get written where” decision.
- Ensure outputs are emitted even when partial results are produced (with explicit failure/skip annotations).

**Must not**
- Compute statistics.
- Re-interpret domain semantics.
- Depend on the source file format.

**Handoff artifact**
- Concrete artifacts on disk: tables with stable schema + run metadata.

> Contract: Reporting never decides *how* a statistic is computed; it only decides *how it is represented and stored*.

---

## Normalization philosophy (without external tools)

Because `scarscape` does not shell out to `bcftools` or similar tools:

- When normalization is required for correctness, it must be implemented in-process using deterministic rules.
- When a normalization step is complex, ambiguous, or not feasible to reproduce faithfully, `scarscape` should:
  1. detect the condition,
  2. produce a clear error or warning with remediation guidance, and
  3. (optionally) provide a mode that assumes pre-normalized inputs.

This preserves the self-contained binary while still being pragmatic about real-world variation in file-formats.

---

## Reference and “bundled resources” philosophy

- Reference-dependent assets (region sets, known normalizations, genome build metadata) are treated as **reference bundles**.
- Bundles can be shipped with the repository and optionally packaged with releases.
- The CLI should allow overriding bundle paths (e.g. `--references-dir`) without requiring code changes.

The reference bundle is the single place to add new region sets and build-specific metadata.

---

## Extension rules (how to add new metrics safely)

To add a new statistic:
1. Extend or reuse an I/O adapter to expose the required normalized event stream(s).
2. Implement the computation in a stats engine module, using existing primitives when possible.
3. Add a table schema (or extend an existing one) in the reporting layer.
4. Wire it into the pipeline orchestration with explicit error/skip behavior.

If adding a metric requires parsing raw VCF fields in the stats layer, the design is wrong. Move parsing/validation up into adapters.

---

## Defaults and UX

- Errors should be actionable and scoped:
  - include sample ID, file path, and the specific violated assumption.
- Partial results are acceptable only when explicitly labeled in outputs.

---

## What “done” looks like

A feature is “done” when:
- it runs deterministically and consumes no more than <2GB of memory,
- it emits schema-stable outputs,
- it fails with actionable messages on malformed inputs,
- it is covered by unit tests at the adapter boundary and engine boundary, and
- it does not introduce cross-layer leakage (e.g. VCF parsing in stats, formatting in engines).
