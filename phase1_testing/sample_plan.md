# Phase 1: MetaGrouper Biological Validation Plan

## Target: 50 samples from 5 environments

### Environment 1: Human Gut (10 samples)
- **Source**: American Gut Project (PRJEB11419)
- **Expected clustering**: By individual/diet
- **Sample size**: 50-200MB each

### Environment 2: Marine Water (10 samples) 
- **Source**: Ocean Sampling Day (OSD)
- **Expected clustering**: By depth/location
- **Sample size**: 50-300MB each

### Environment 3: Soil (10 samples)
- **Source**: EMP soil studies
- **Expected clustering**: By pH/chemistry
- **Sample size**: 100-400MB each

### Environment 4: Freshwater (10 samples)
- **Source**: Lake microbiome studies
- **Expected clustering**: By lake/depth
- **Sample size**: 50-200MB each

### Environment 5: Built Environment (10 samples)
- **Source**: Built environment microbiome
- **Expected clustering**: By location type
- **Sample size**: 30-150MB each

## Total Expected Size: 5-15GB

## Download Strategy
1. Start with smallest samples from each environment
2. Use `fastq-dump --split-files` for SRA format
3. Test MetaGrouper on each environment separately first
4. Combine all environments for final validation

## Success Criteria
- Clear clustering by environment type
- Sensible assembly recommendations within environments
- All phases (1-4) complete without errors
- Processing time < 30 minutes for 50 samples

## Backup Plan
If any environment fails to download or cluster properly:
- Replace with alternative samples from same environment type
- Minimum 30 samples from 3 environments still validates the approach