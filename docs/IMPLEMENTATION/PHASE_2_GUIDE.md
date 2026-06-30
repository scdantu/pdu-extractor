# Phase 2: Feature Extraction Implementation Guide

## Overview
Implement the feature extraction pipeline as modular, reusable classes that support multiple encoding strategies (20-AA, functional grouping, conservation metrics).

## Files to Create

### 1. `pdusearch/features/encodings.py` (Core)

```python
from enum import Enum
from abc import ABC, abstractmethod
from typing import Tuple, List
import numpy as np

class Encoding(Enum):
    """Feature encoding strategies."""
    AA_20 = "aa_20"                              # Current: 20 AA types
    FUNCTIONAL_5 = "functional_5"                # Task #1: 5 functional classes
    AA_20_WITH_CONSERVATION = "aa_20_conservation"  # Task #2: 20 AA + entropy

class EncodingStrategy(ABC):
    """Abstract base for encoding strategies."""
    
    @abstractmethod
    def get_n_features(self, n_distance_bins: int) -> int:
        """Return feature dimension."""
        pass
    
    @abstractmethod
    def encode_neighbors(self, neighbors: List) -> np.ndarray:
        """Convert neighbor data to feature vector."""
        pass

class AA20Encoding(EncodingStrategy):
    """Standard 20-amino-acid encoding."""
    
    def get_n_features(self, n_distance_bins: int) -> int:
        return 20 * 3 * n_distance_bins  # 20 AAs × 3 SS × distance bins
    
    def encode_neighbors(self, neighbors):
        # Build 900D vector from neighbor list
        # neighbors = [(residue_type, ss, distance), ...]
        pass

class Functional5Encoding(EncodingStrategy):
    """Functional class encoding (hydrophobic, polar, charged+/-, special)."""
    
    def get_n_features(self, n_distance_bins: int) -> int:
        return 5 * 3 * n_distance_bins  # 5 functional classes × 3 SS × distance bins
    
    def encode_neighbors(self, neighbors):
        # Group residues by FunctionalClass, build 225D vector
        pass

class AA20ConservationEncoding(EncodingStrategy):
    """20-AA encoding with Shannon entropy per position."""
    
    def get_n_features(self, n_distance_bins: int) -> int:
        return 20 * 3 * n_distance_bins + n_distance_bins  # +entropy per bin
    
    def encode_neighbors(self, neighbors, conservation_scores=None):
        # Build 900D vector + add conservation metrics
        pass

# Factory function
def get_encoding_strategy(encoding: Encoding) -> EncodingStrategy:
    """Get strategy for encoding type."""
    strategies = {
        Encoding.AA_20: AA20Encoding(),
        Encoding.FUNCTIONAL_5: Functional5Encoding(),
        Encoding.AA_20_WITH_CONSERVATION: AA20ConservationEncoding(),
    }
    return strategies[encoding]
```

**Key design**: Each encoding is independent, easy to test, easy to extend to new types.

### 2. `pdusearch/features/alignment.py` (PyCoM Integration)

```python
from typing import Dict, List, Optional
import numpy as np

class AlignmentProvider(ABC):
    """Abstract alignment data provider."""
    
    @abstractmethod
    def get_alignment(self, pdb_id: str, chain_id: str) -> np.ndarray:
        """Get alignment matrix (n_positions × n_residues_observed)."""
        pass
    
    @abstractmethod
    def get_conservation_score(self, pdb_id: str, chain_id: str, position: int) -> float:
        """Get Shannon entropy or conservation score."""
        pass

class PyComAlignmentProvider(AlignmentProvider):
    """PyCoM-based alignment provider."""
    
    def __init__(self, pycom_db: str):
        # Initialize PyCoM database connection
        pass
    
    def get_alignment(self, pdb_id: str, chain_id: str):
        # Query PyCoM HHBLITS alignments
        pass
    
    def get_conservation_score(self, pdb_id: str, chain_id: str, position: int):
        # Compute Shannon entropy from alignment
        pass

def compute_shannon_entropy(alignment_column: np.ndarray) -> float:
    """Compute Shannon entropy of alignment column."""
    counts = np.bincount(alignment_column)
    probs = counts / len(alignment_column)
    entropy = -np.sum(probs[probs > 0] * np.log2(probs[probs > 0]))
    return entropy
```

### 3. `pdusearch/features/extractor.py` (Main Class)

```python
from pdusearch.config import Config, FeatureConfig
from pdusearch.io import PDUDatabase, BatchPDUFetcher
from pdusearch.utils import AA_ORDER, SS_ORDER, AA_TO_FUNCTIONAL
from .encodings import Encoding, get_encoding_strategy
from .alignment import AlignmentProvider, PyComAlignmentProvider
import numpy as np

class FeatureExtractor:
    """Extract PDU features with configurable encoding."""
    
    def __init__(
        self,
        config: Config,
        encoding: Encoding = Encoding.AA_20,
        alignment_provider: Optional[AlignmentProvider] = None,
    ):
        """Initialize extractor.
        
        Args:
            config: Config instance with db path, radius, etc.
            encoding: Encoding strategy to use
            alignment_provider: Optional alignment provider for conservation metrics
        """
        self.config = config
        self.encoding = encoding
        self.alignment_provider = alignment_provider
        
        self.db = PDUDatabase(config.db)
        self.strategy = get_encoding_strategy(encoding)
        
        self.logger = get_logger(__name__)
    
    def extract_for_aa(self, aa: str, max_distance: Optional[float] = None):
        """Extract features for single amino acid.
        
        Args:
            aa: Amino acid one-letter code (L, A, G, etc.)
            max_distance: Override config radius if specified
        
        Returns:
            (X, pdu_ids) where X is (n_pdus, n_features) and pdu_ids is (n_pdus,)
        """
        distance = max_distance or self.config.radius
        
        self.logger.info(f"Extracting features for {aa} at {distance}Å...")
        
        with self.db as db:
            pdu_ids = db.get_pdu_ids_for_aa(aa)
            neighbors_data = db.get_pdu_neighbors(aa, max_distance=distance)
        
        # Build feature matrix
        n_features = self.strategy.get_n_features(self.config.n_distance_bins)
        X = np.zeros((len(pdu_ids), n_features), dtype=np.float32)
        
        # For each PDU, encode its neighbors
        for idx, pdu_id in enumerate(pdu_ids):
            neighbors = [n for n in neighbors_data if n[0] == pdu_id]
            X[idx] = self.strategy.encode_neighbors(neighbors)
        
        self.logger.info(f"  → Extracted {len(pdu_ids):,} PDUs × {n_features} features")
        
        return X, np.array(pdu_ids)
    
    def extract_batch(self, aa_list: List[str]):
        """Extract features for multiple amino acids.
        
        Returns:
            Dict mapping aa -> (X, pdu_ids)
        """
        results = {}
        for aa in aa_list:
            results[aa] = self.extract_for_aa(aa)
        return results
    
    def save_features(self, X: np.ndarray, pdu_ids: np.ndarray, output_path: str):
        """Save features to NPZ file."""
        np.savez(
            output_path,
            X=X,
            pdu_ids=pdu_ids,
            encoding=self.encoding.value,
            radius=self.config.radius,
        )
        self.logger.info(f"Saved to {output_path}")
```

## Implementation Order

1. **`encodings.py`** first (no dependencies)
   - Define Encoding enum
   - Implement AA20Encoding (simplest, copy from export_pdu_features.py)
   - Implement Functional5Encoding (group AAs by FunctionalClass)
   - Implement AA20ConservationEncoding (placeholder, Task #2)

2. **`alignment.py`** second (independent)
   - Define AlignmentProvider ABC
   - Implement PyComAlignmentProvider (Task #6)
   - Implement conservation computation

3. **`extractor.py`** last (builds on 1 & 2)
   - FeatureExtractor class
   - Integration with PDUDatabase
   - Batch operations

## Testing Strategy

### Test 1: Encoding Strategies
```python
from pdusearch.features.encodings import AA20Encoding, Functional5Encoding

encoding_aa20 = AA20Encoding()
assert encoding_aa20.get_n_features(15) == 900

encoding_func = Functional5Encoding()
assert encoding_func.get_n_features(15) == 225
```

### Test 2: Feature Extraction
```python
from pdusearch.features import FeatureExtractor, Encoding
from pdusearch.config import Config

config = Config(db="per_aa_sqlite/pdus_L.sqlite", aa="L")
extractor = FeatureExtractor(config, encoding=Encoding.AA_20)
X, pdu_ids = extractor.extract_for_aa("L")

assert X.shape == (n_pdus_for_L, 900)
assert len(pdu_ids) == n_pdus_for_L
```

### Test 3: Multiple Encodings
```python
X_aa20, ids = extractor_aa20.extract_for_aa("L")
X_func, _ = extractor_func.extract_for_aa("L")

assert X_aa20.shape[0] == X_func.shape[0]  # Same PDUs
assert X_aa20.shape[1] == 900  # 20 AA × 3 SS × 15 bins
assert X_func.shape[1] == 225  # 5 func × 3 SS × 15 bins
```

## Integration with Existing Code

### Migration Path
1. Keep existing `bin/extract/export_pdu_features.py` as reference
2. New `bin/extract_features.py` calls `FeatureExtractor` from library
3. Eventually deprecate old script

### CLI Wrapper Example
```python
#!/usr/bin/env python3
"""Extract PDU features (new library-based version)."""

import argparse
from pdusearch.features import FeatureExtractor, Encoding
from pdusearch.config import Config
from pdusearch.logging_utils import configure_logging, add_logging_args

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default="per_aa_sqlite/pdus_L.sqlite")
    parser.add_argument("--aa", required=True)
    parser.add_argument("--encoding", default="aa_20", choices=["aa_20", "functional_5", "aa_20_conservation"])
    parser.add_argument("--radius", type=float, default=15.0)
    parser.add_argument("--out-dir", default="analysis/features")
    add_logging_args(parser)
    args = parser.parse_args()
    
    configure_logging(args.log_file, args.log_level)
    
    config = Config(db=args.db, aa=args.aa, radius=args.radius, out_dir=args.out_dir)
    extractor = FeatureExtractor(config, encoding=Encoding[args.encoding.upper()])
    X, ids = extractor.extract_for_aa(args.aa)
    
    output_file = f"{args.out_dir}/pdu_features_{args.aa}_{args.encoding}.npz"
    extractor.save_features(X, ids, output_file)

if __name__ == "__main__":
    main()
```

## Dependencies

All imports are already available:
- ✓ `numpy` - Feature arrays
- ✓ `sqlite3` - Database (via PDUDatabase)
- ✓ `dataclasses` - Config
- ✓ `enum` - Encoding types
- ✓ `logging` - Logging
- Later: `pycom` library (Task #6), `scipy` for statistics (Task #2)

---

## Deliverables for Phase 2

1. ✓ `pdusearch/features/encodings.py` - 3 encoding strategies
2. ✓ `pdusearch/features/alignment.py` - PyCoM integration framework
3. ✓ `pdusearch/features/extractor.py` - Main FeatureExtractor class
4. ✓ `pdusearch/features/__init__.py` - Export public API
5. ✓ Tests for each module
6. ✓ `bin/extract_features.py` - CLI wrapper

Once Phase 2 is done:
- Task #1 (functional categories) - ✓ Encoding.FUNCTIONAL_5 ready
- Task #2 (conservation metrics) - ✓ Framework ready, PyCoM integration pending Task #6
- Task #6 (PyCoM alignment data) - Can proceed immediately with existing framework
