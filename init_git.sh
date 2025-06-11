#!/bin/bash

# MetaGrouper Git Repository Initialization Script
# This script initializes the git repository and prepares for first commit

echo "🚀 Initializing MetaGrouper Git Repository"
echo "=========================================="

# Check if we're already in a git repository
if [ -d ".git" ]; then
    echo "⚠️  Git repository already exists. Skipping initialization."
else
    echo "📁 Initializing git repository..."
    git init
fi

# Add all files to staging (respecting .gitignore)
echo "📋 Adding files to staging area..."
git add .

# Show what will be committed
echo "📊 Files to be committed:"
git status --short

# Create initial commit
echo "💾 Creating initial commit..."
git commit -m "Initial commit: Complete MetaGrouper implementation

- Phase 1: K-mer profiling and sample similarity analysis
- Phase 2: Metadata variable analysis with PERMANOVA  
- Phase 3: Assembly strategy recommendations
- Phase 4: CLI interface and comprehensive documentation
- Complete test suite and examples
- Production-ready package structure"

echo ""
echo "✅ Git repository initialized successfully!"
echo ""
echo "📋 Next steps:"
echo "1. Add your GitHub repository as remote:"
echo "   git remote add origin git@github.com:megjohnson1999/metaGrouper.git"
echo ""
echo "2. Push to GitHub:"
echo "   git branch -M main"
echo "   git push -u origin main"
echo ""
echo "3. Verify everything looks good on GitHub"
echo ""
echo "🎉 MetaGrouper is ready for GitHub!"