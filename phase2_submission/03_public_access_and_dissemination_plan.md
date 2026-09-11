# OligoTox Open Data Challenge — Public Access and Dissemination Plan (PADP)

**Submitter:** Anote, Inc. | **Point of Contact:** Natan Vidra
**Submission target:** Single PDF, ≤5 pages, 8.5"×11", ≥1" margins, ≥11pt Arial, line spacing ≥1.0

> Per the Rules section of the Challenge Announcement, winning this prize requires agreeing to
> abide by the terms of this PADP, and NIH intends to publicly post it. This document describes
> *commitments*, not experimental results, so it does not depend on the dataset's final size or
> composition — only the repository names in Section 2 should be double-checked against what was
> actually used once the dataset is deposited.

## 1. Overview

Anote, Inc. commits to releasing all data generated under this Challenge — compound registry,
computed chemical/structural features, raw and normalized experimental results, and derived
summary statistics — as an open, FAIR-compliant dataset with no proprietary restrictions, no
embargo beyond a short rolling-release schedule, and no registration barrier to download.

## 2. Dissemination of the Solution and the Knowledge Needed to Use It

*(Addresses: "how winners will disseminate information about the solution and make the
solution... available under non-exclusive licenses for research purposes.")*

- **Primary data repository:** Zenodo (a persistent DOI is assigned automatically at deposit)
  for compound registry, computed features, assay results, and derived statistics, in
  CSV/Parquet format.
- **Domain-specific repositories:** raw transcriptomics to GEO (NCBI), raw proteomics to PRIDE
  Archive (EBI), each in the community-standard raw format (FASTQ, mzML respectively).
- **Data dictionary and schema:** published alongside the dataset (see
  `dataset/data_dictionary_and_schema.md`) documenting every field, unit, and controlled
  vocabulary, so the dataset is usable without needing to consult Anote directly.
- **Code and tooling:** ingestion, QC, feature-extraction, and modeling code released on GitHub
  (this repository) under Apache 2.0.
- **Non-exclusive licensing:** all data released under **Creative Commons Attribution 4.0 (CC
  BY 4.0)** — permits any use, including commercial, with attribution only, for any researcher
  without a bilateral agreement with Anote. All code released under **Apache 2.0**.
- **Findability:** persistent DOIs, metadata registered with DataCite/re3data, and — where
  applicable — indexing in domain registries (e.g., FAIRsharing.org) so the dataset is
  discoverable by researchers who have never interacted with Anote directly.
- **Direct outreach:** notification to relevant scientific communities (e.g., the
  Oligonucleotide Therapeutics Society) and a public data-descriptor writeup accompanying the
  release, so uptake does not depend on individual requests to Anote.

## 3. Continuity Plan if Anote Is Unable to Maximize Public Access

*(Addresses: "how the winner would allow others to utilize the solution in the event the winner
is unable themselves to maximize public access... specific licensing schemes and scenarios when
such schemes would be employed.")*

Should Anote become unable to maintain active distribution (e.g., company dissolution,
discontinued hosting, loss of maintainers), the following scenarios and licensing terms apply:

- **Scenario A — Hosting lapses, data already released:** Because the primary release sits on
  third-party, long-term-preservation repositories (Zenodo/EMBL-EBI commit to multi-decade
  archival; GEO/PRIDE are NIH/NCBI- and EBI-maintained), the CC BY 4.0 / Apache 2.0 licenses
  already granted remain in force and fully usable by anyone, independent of Anote's continued
  operation. No further action is required from Anote for the data to stay accessible.
- **Scenario B — Data generated but not yet released:** If any portion of the dataset exists
  but has not yet been deposited at the time Anote becomes unable to continue, Anote grants, in
  advance, an irrevocable, worldwide, royalty-free, non-exclusive license under CC BY 4.0 (data)
  / Apache 2.0 (code) to any successor maintainer, research consortium, or archival body (e.g.,
  a repository operated by NCATS, NCBI, or a designated academic partner) to receive, host, and
  redistribute the unreleased material on Anote's behalf, so publication is not blocked on
  Anote's continued existence.
- **Scenario C — Maintenance capacity lapses but Anote still exists:** Anote will transfer
  repository/portal administrative access to a mutually agreed academic or nonprofit steward
  (e.g., a university library data-curation group) under the same open licenses, rather than
  allow the resource to go stale.

## 4. Fallback: Government-Enabled Access if Anote Fails and Does Not License Others

*(Addresses: "how the winner will permit the U.S. government to allow interested parties to
utilize the solution if the winner themselves fails to utilize the solution and does not permit
others to utilize the solution under reasonable terms.")*

Independent of, and in addition to, the nonexclusive license already granted to NIH under the
Challenge Rules (to reproduce, publish, post, link to, share, and display the submission,
including this PADP, and to practice or have practiced the underlying solution), Anote commits
that:

- If Anote fails to release the dataset/code as committed above, or refuses to permit third
  parties to use it under reasonable (i.e., CC BY 4.0 / Apache 2.0-equivalent) terms, Anote's
  existing license grant to the federal government is sufficient for NIH/NCATS to exercise its
  rights to make the solution (data, schema, and any deposited code) available to interested
  researchers directly — including by depositing Anote-provided raw data files with NIH/NCATS
  at time of submission (see `dataset/README_raw_data_access.md`) as a standing backup copy, so
  the government is never solely dependent on Anote's own infrastructure to exercise this right.
- Anote will not impose additional restrictions, paywalls, or exclusive licenses on this data at
  any point after award that would conflict with the open terms described in Section 2.

## 5. Licensing Summary

| Asset | License |
|---|---|
| Experimental data (compound registry, results, features) | CC BY 4.0 |
| Data dictionary / schema | CC BY 4.0 |
| Processing/ingestion/QC/model code | Apache 2.0 |
| Analysis notebooks | Apache 2.0 |

## 6. Long-Term Sustainability

- Long-term archival hosts (Zenodo, GEO, PRIDE) chosen specifically for multi-decade
  preservation commitments independent of any single organization's continued operation.
- Anote commits to maintaining active repositories/tooling for a minimum of 3 years post-award,
  with the governance/maintenance transfer plan (Section 3, Scenario C) as a backstop beyond
  that, on top of the Zenodo/GEO/PRIDE archival commitments in Section 3, Scenario A, which do
  not depend on Anote's continued operation at all.

---

*Prepared by Natan Vidra, Anote, Inc. Contact: nvidra@anote.ai*
