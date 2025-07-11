# MetaGrouper Enhancements

## New Features Added

### 1. Phase Selection Support (`metagrouper_phases.py`)

**New command-line options:**
- `--phases 1 2 3` - Run only specific phases
- `--skip-phases 1` - Skip specific phases  
- `--load-from path/` - Load Phase 1 results from previous run

**Benefits:**
- Skip expensive k-mer computation when reanalyzing with different metadata
- Debug specific phases without running everything
- Try different parameters for Phase 2/3 without recomputing profiles
- Save time and computational resources

**Usage Examples:**
```bash
# Run only k-mer profiling (Phase 1)
python metagrouper_phases.py /path/to/fastq/files -o profiles/ --phases 1

# Rerun metadata analysis with different settings
python metagrouper_phases.py dummy_path -m new_metadata.csv -o reanalysis/ \
  --load-from profiles/ --phases 2

# Run phases 2 and 3 using existing profiles
python metagrouper_phases.py dummy_path -m metadata.csv -o analysis/ \
  --load-from profiles/ --skip-phases 1
```

### 2. Metadata Loading Bug Fix

**Problem:** MetaGrouper failed to analyze any metadata variables, incorrectly classifying all variables as "essentially constant values."

**Root Cause:** Type mismatch between FASTQ sample names (strings like `'10571'`) and metadata database_ID values (integers like `10571`). This caused pandas reindexing to fail, creating all-NaN rows.

**Fix:** Enhanced `load_metadata()` method in `MetadataAnalyzer` class to:
1. Detect metadata index data type
2. Convert sample names to match metadata type (strings to integers)
3. Fall back to converting metadata index to strings if conversion fails

**Files Modified:**
- `phases/metadata_analyzer.py` - Lines 471-491

**Result:** PERMANOVA analysis now works correctly with integer database IDs and string sample names.

## Testing

✅ **Phase Selection Logic:** Verified all argument combinations work correctly
✅ **Metadata Loading Fix:** Tested with sample data matching the original issue
✅ **Backward Compatibility:** Original functionality preserved
✅ **Help Documentation:** Comprehensive usage examples included

## Deployment

To apply these changes:

1. **Direct commit:** Add and commit the new files
2. **Patch file:** Apply `metagrouper_enhancements.patch`

## Impact

For the original 1150 gut microbiome sample analysis:
- **Previous:** 0/56 metadata variables analyzed (all excluded as "constant")
- **After fix:** All valid variables (IBD, sex, diagnosis, etc.) correctly analyzed
- **Time savings:** Can rerun Phase 2 in minutes instead of hours by reusing k-mer profiles

## Files Added/Modified

- **NEW:** `metagrouper_phases.py` - Enhanced version with phase selection
- **MODIFIED:** `phases/metadata_analyzer.py` - Bug fix for type mismatch
- **NEW:** `CHANGES.md` - This documentation
- **NEW:** `metagrouper_enhancements.patch` - Patch file for easy deployment