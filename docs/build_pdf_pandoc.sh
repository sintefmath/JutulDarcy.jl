#!/bin/bash
# Build PDF documentation for JutulDarcy.jl using Pandoc
# This script collects all markdown documentation and compiles it to PDF

set -e  # Exit on error

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}Building JutulDarcy.jl PDF Documentation${NC}"
echo "============================================="

# Check if pandoc is installed
if ! command -v pandoc &> /dev/null; then
    echo -e "${RED}Error: pandoc is not installed${NC}"
    echo "Please install pandoc from https://pandoc.org/installing.html"
    exit 1
fi

# Check if xelatex is installed (required by pandoc for PDF output with Unicode)
if ! command -v xelatex &> /dev/null; then
    echo -e "${YELLOW}Warning: xelatex is not installed${NC}"
    echo "Pandoc requires a LaTeX engine for PDF output."
    echo "Please install a LaTeX distribution:"
    echo "  - Linux: sudo apt-get install texlive-xetex"
    echo "  - macOS: brew install --cask mactex"
    echo "  - Windows: Install MiKTeX or TeX Live"
    exit 1
fi

echo -e "${GREEN}Prerequisites check passed${NC}"

# Create temporary directory for processing
TMP_DIR=$(mktemp -d)
echo "Using temporary directory: $TMP_DIR"

# Copy the PDF-friendly index
if [ -f "src/index_pdf.md" ]; then
    cp "src/index_pdf.md" "$TMP_DIR/00_index.md"
else
    echo -e "${YELLOW}Warning: src/index_pdf.md not found, using regular index${NC}"
    # Remove Vitepress frontmatter from index.md
    sed '/^````@raw html$/,/^````$/d' "src/index.md" > "$TMP_DIR/00_index.md"
fi

# Function to clean markdown files (remove Documenter-specific syntax)
clean_markdown() {
    local input="$1"
    local output="$2"
    
    # Remove @ref, @docs, @example blocks and other Documenter-specific syntax
    sed -e 's/@ref[[:space:]]*[^)]*//g' \
        -e 's/@docs//g' \
        -e '/^```@example/,/^```$/d' \
        -e '/^```@raw html$/,/^```$/d' \
        -e 's/@raw html//g' \
        "$input" > "$output"
}

# Generate example overview first
echo "Generating example overview..."
julia generate_example_overview.jl

# Collect pages in the correct order following make.jl structure
echo "Collecting documentation pages..."
COUNTER=1

# 1. Introduction section (from index_pdf.md)
if [ -f "src/index_pdf.md" ]; then
    printf -v padded "%02d" $COUNTER
    clean_markdown "src/index_pdf.md" "$TMP_DIR/${padded}_index.md"
    echo "  Added: Introduction"
    COUNTER=$((COUNTER + 1))
fi

# 2. Getting started
if [ -f "src/man/intro.md" ]; then
    printf -v padded "%02d" $COUNTER
    clean_markdown "src/man/intro.md" "$TMP_DIR/${padded}_intro.md"
    echo "  Added: Getting started"
    COUNTER=$((COUNTER + 1))
fi

# 3. First example
if [ -f "src/man/first_ex.md" ]; then
    printf -v padded "%02d" $COUNTER
    clean_markdown "src/man/first_ex.md" "$TMP_DIR/${padded}_first_ex.md"
    echo "  Added: Your first simulation"
    COUNTER=$((COUNTER + 1))
fi

# 4. FAQ
if [ -f "src/extras/faq.md" ]; then
    printf -v padded "%02d" $COUNTER
    clean_markdown "src/extras/faq.md" "$TMP_DIR/${padded}_faq.md"
    echo "  Added: FAQ"
    COUNTER=$((COUNTER + 1))
fi

# 5. Fundamentals section
echo "  Section: Fundamentals"
for section in highlevel basics/input_files basics/systems basics/solution; do
    if [ -f "src/man/$section.md" ]; then
        printf -v padded "%02d" $COUNTER
        clean_markdown "src/man/$section.md" "$TMP_DIR/${padded}_$(basename $section).md"
        echo "    Added: man/$section.md"
        COUNTER=$((COUNTER + 1))
    fi
done

# 6. Detailed API section
echo "  Section: Detailed API"
for section in basics/forces basics/wells basics/primary basics/secondary basics/parameters basics/plotting basics/utilities; do
    if [ -f "src/man/$section.md" ]; then
        printf -v padded "%02d" $COUNTER
        clean_markdown "src/man/$section.md" "$TMP_DIR/${padded}_$(basename $section).md"
        echo "    Added: man/$section.md"
        COUNTER=$((COUNTER + 1))
    fi
done

# 7. Parallelism and compilation section
echo "  Section: Parallelism and compilation"
for section in advanced/mpi advanced/gpu advanced/compiled; do
    if [ -f "src/man/$section.md" ]; then
        printf -v padded "%02d" $COUNTER
        clean_markdown "src/man/$section.md" "$TMP_DIR/${padded}_$(basename $section).md"
        echo "    Added: man/$section.md"
        COUNTER=$((COUNTER + 1))
    fi
done

# 8. References section
echo "  Section: References"
if [ -f "src/man/basics/package.md" ]; then
    printf -v padded "%02d" $COUNTER
    clean_markdown "src/man/basics/package.md" "$TMP_DIR/${padded}_package.md"
    echo "    Added: man/basics/package.md"
    COUNTER=$((COUNTER + 1))
fi

if [ -f "src/extras/paper_list.md" ]; then
    printf -v padded "%02d" $COUNTER
    clean_markdown "src/extras/paper_list.md" "$TMP_DIR/${padded}_paper_list.md"
    echo "    Added: extras/paper_list.md"
    COUNTER=$((COUNTER + 1))
fi

if [ -f "src/ref/jutul.md" ]; then
    printf -v padded "%02d" $COUNTER
    clean_markdown "src/ref/jutul.md" "$TMP_DIR/${padded}_jutul_ref.md"
    echo "    Added: ref/jutul.md"
    COUNTER=$((COUNTER + 1))
fi

if [ -f "src/extras/refs.md" ]; then
    printf -v padded "%02d" $COUNTER
    clean_markdown "src/extras/refs.md" "$TMP_DIR/${padded}_refs.md"
    echo "    Added: extras/refs.md"
    COUNTER=$((COUNTER + 1))
fi

# 9. Examples section
echo "  Section: Examples"
if [ -f "src/examples/overview/example_overview.md" ]; then
    printf -v padded "%02d" $COUNTER
    clean_markdown "src/examples/overview/example_overview.md" "$TMP_DIR/${padded}_example_overview.md"
    echo "    Added: examples/overview/example_overview.md"
    COUNTER=$((COUNTER + 1))
fi

# 10. Validation section
echo "  Section: Validation"
if [ -f "src/man/validation.md" ]; then
    printf -v padded "%02d" $COUNTER
    clean_markdown "src/man/validation.md" "$TMP_DIR/${padded}_validation.md"
    echo "    Added: man/validation.md"
    COUNTER=$((COUNTER + 1))
fi

# Combine all markdown files
echo "Combining markdown files..."
COMBINED_MD="$TMP_DIR/combined.md"
cat $TMP_DIR/*.md > "$COMBINED_MD"

# Build PDF with pandoc
OUTPUT_PDF="JutulDarcy_Documentation.pdf"
echo -e "${GREEN}Building PDF with Pandoc...${NC}"
echo "This may take a few minutes..."

pandoc "$COMBINED_MD" \
    -o "$OUTPUT_PDF" \
    --pdf-engine=xelatex \
    --toc \
    --toc-depth=3 \
    --number-sections \
    --highlight-style=tango \
    --variable=geometry:margin=1in \
    --variable=documentclass:report \
    --variable=fontsize:11pt \
    --variable=papersize:letter \
    --metadata title="JutulDarcy.jl Documentation" \
    --metadata author="Olav Møyner and contributors" \
    --metadata date="$(date +%Y-%m-%d)" \
    2>&1 | grep -v "Missing character" || true

# Clean up
echo "Cleaning up temporary files..."
rm -rf "$TMP_DIR"

if [ -f "$OUTPUT_PDF" ]; then
    echo -e "${GREEN}=============================================${NC}"
    echo -e "${GREEN}Success!${NC}"
    echo "PDF documentation created: $OUTPUT_PDF"
    echo "File size: $(du -h $OUTPUT_PDF | cut -f1)"
else
    echo -e "${RED}Error: PDF generation failed${NC}"
    exit 1
fi
