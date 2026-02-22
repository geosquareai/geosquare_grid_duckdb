#!/bin/bash
# Script to update DuckDB extension to a specific version or latest version

set -e

REPO_DIR="/mnt/data/REPO/geosquare_lib/duckdb"
TARGET_VERSION=$1

if [ -z "$TARGET_VERSION" ]; then
    # Get latest tag from DuckDB repo
    TARGET_VERSION=$(git ls-remote --tags https://github.com/duckdb/duckdb.git | grep -o 'refs/tags/v[0-9]*\.[0-9]*\.[0-9]*' | cut -d/ -f3 | sort -V | tail -n 1)
    echo "No version specified. Latest version found: $TARGET_VERSION"
fi

echo "Updating to DuckDB version $TARGET_VERSION..."

# Update submodules
cd "$REPO_DIR"
git submodule update --init --recursive
cd duckdb
git fetch --tags
git checkout "$TARGET_VERSION"
cd ..
cd extension-ci-tools
git fetch --tags
# Note: extension-ci-tools might not always have the exact same tag, usually they match
if git rev-parse "$TARGET_VERSION" >/dev/null 2>&1; then
    git checkout "$TARGET_VERSION"
else
    echo "Warning: $TARGET_VERSION not found in extension-ci-tools, staying on current branch."
fi
cd ..

# Update CI workflow
WORKFLOW_FILE=".github/workflows/MainDistributionPipeline.yml"
if [ -f "$WORKFLOW_FILE" ]; then
    # Replace old version with new version in the duckdb-stable-build job
    # This regex is a bit simple, might need adjustment depending on the file structure
    sed -i "s/duckdb_version: v[0-9]*\.[0-9]*\.[0-9]*/duckdb_version: $TARGET_VERSION/g" "$WORKFLOW_FILE"
    sed -i "s/ci_tools_version: v[0-9]*\.[0-9]*\.[0-9]*/ci_tools_version: $TARGET_VERSION/g" "$WORKFLOW_FILE"
    sed -i "s/@v[0-9]*\.[0-9]*\.[0-9]*/@$TARGET_VERSION/g" "$WORKFLOW_FILE"
    echo "Updated $WORKFLOW_FILE"
fi

echo "Successfully updated to $TARGET_VERSION. Please commit and push the changes."
