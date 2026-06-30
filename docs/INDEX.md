# Documentation Index

## Getting Started
- **[README](../README.md)** — Project overview and quick start
- **[Quick Start](GUIDES/QUICK_START.md)** — Installation and first run
- **[Testing Guide](TESTING/COMPRESSION_TEST_GUIDE.md)** — Run compression hypothesis test

## Architecture & Implementation
- **[Architecture Overview](IMPLEMENTATION/REFACTORING_COMPLETE.md)** — Library design and structure
- **[Refactoring Plan](IMPLEMENTATION/REFACTOR_PLAN.md)** — Original planning document

### Implementation Phases
- **[Phase 1: Core Infrastructure](IMPLEMENTATION/PHASE_1_SUMMARY.md)** — Config, logging, I/O
- **[Phase 2: Feature Extraction](IMPLEMENTATION/PHASE_2_SUMMARY.md)** — 3 encoding strategies
  - [Phase 2 Guide](IMPLEMENTATION/PHASE_2_GUIDE.md) — Detailed implementation guide
  - [Phase 2 Checklist](IMPLEMENTATION/PHASE_2_CHECKLIST.md) — Component verification
- **[Phase 3: Model Training](IMPLEMENTATION/PHASE_3_SUMMARY.md)** — Autoencoders

## Testing

- **[Compression Test Guide](TESTING/COMPRESSION_TEST_GUIDE.md)** — Complete testing documentation
- **[Implementation Status](IMPLEMENTATION/README_IMPLEMENTATION.md)** — Current phase status

## API Reference

Auto-generated documentation from code:
- [pdusearch.features](../pdusearch/features/__init__.py) — Feature extraction
- [pdusearch.models](../pdusearch/models/__init__.py) — Model training
- [pdusearch.clustering](../pdusearch/clustering/__init__.py) — Clustering
- [pdusearch.config](../pdusearch/config.py) — Configuration classes

## Quick Reference

| Need | Document |
|------|-----------|
| Get started quickly | [Quick Start](GUIDES/QUICK_START.md) |
| Understand architecture | [Architecture Overview](IMPLEMENTATION/REFACTORING_COMPLETE.md) |
| Run compression test | [Testing Guide](TESTING/COMPRESSION_TEST_GUIDE.md) |
| Feature extraction API | [Phase 2 Summary](IMPLEMENTATION/PHASE_2_SUMMARY.md) |
| Model training API | [Phase 3 Summary](IMPLEMENTATION/PHASE_3_SUMMARY.md) |
| Implementation details | [Implementation Phases](IMPLEMENTATION/) |

## File Organization

```
docs/
├── INDEX.md                    # This file
├── IMPLEMENTATION/
│   ├── REFACTOR_PLAN.md       # Original planning
│   ├── REFACTORING_COMPLETE.md # Architecture & design
│   ├── PHASE_1_SUMMARY.md     # Core infrastructure
│   ├── PHASE_2_GUIDE.md       # Feature extraction guide
│   ├── PHASE_2_SUMMARY.md     # Feature extraction summary
│   ├── PHASE_2_CHECKLIST.md   # Verification checklist
│   ├── PHASE_3_SUMMARY.md     # Model training
│   └── README_IMPLEMENTATION.md # Status & next steps
├── TESTING/
│   └── COMPRESSION_TEST_GUIDE.md # Testing documentation
└── GUIDES/
    └── QUICK_START.md          # Setup & first run
```

---

**Quick Navigation**:
- [Back to README](../README.md)
- [Start Testing](TESTING/COMPRESSION_TEST_GUIDE.md)
- [View Architecture](IMPLEMENTATION/REFACTORING_COMPLETE.md)
