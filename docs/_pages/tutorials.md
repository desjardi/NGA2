---
title: "Tutorials"
permalink: /tutorials/
layout: single
sidebar:
  nav: "tutorials"
---

This page is the tutorial hub for NGA2, starting from “it builds” and moving toward more realistic workflows.

Use the sidebar to navigate through the sequence, or start here:

- [Installing dependencies]({{ '/tutorials/install-deps/' | relative_url }})
- [Installing NGA2]({{ '/tutorials/install-nga2/' | relative_url }})
- [A first run]({{ '/tutorials/first-run/' | relative_url }})
- [A second run]({{ '/tutorials/second-run/' | relative_url }})

## Before you start

- Follow the installation steps:
  - [Installing dependencies]({{ '/tutorials/install-deps/' | relative_url }})
  - [Installing NGA2]({{ '/tutorials/install-nga2/' | relative_url }})
- Pick one example case from the repository’s `examples/` directory.

## Tutorials

### Tutorial 1 — Build and run an example

Goal: compile NGA2 and run a small included example end-to-end.

What you should learn:

- Where example inputs live
- How to launch a run (MPI vs. serial)
- Where outputs/restarts go

### Tutorial 2 — Change one physical parameter

Goal: modify a single parameter in an example input deck and confirm the run changes in an expected way.

What you should learn:

- How to identify key runtime parameters
- How to compare two runs (plots, norms, logs)

### Tutorial 3 — Change resolution and check cost/accuracy

Goal: run the same setup at two resolutions and compare runtime and solution differences.

What you should learn:

- How to change resolution safely
- What to watch for in performance scaling and stability

## Next

- Browse results and media: [Gallery]({{ '/gallery/' | relative_url }})
- Background, people/funding, related repos: [About]({{ '/about/' | relative_url }})
