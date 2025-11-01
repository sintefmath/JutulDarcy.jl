# Building PDF Documentation for JutulDarcy.jl

This directory contains scripts to build both the web-based documentation and a compiled PDF version.

## Building the PDF Documentation

### Prerequisites

1. **Julia** (already installed if you're using this package)
2. **LaTeX Distribution** (required for PDF compilation):
   - **Linux**: Install texlive
     ```bash
     sudo apt-get install texlive-full
     ```
   - **macOS**: Install MacTeX
     ```bash
     brew install --cask mactex
     ```
   - **Windows**: Install MiKTeX or TeX Live
   - **Using Docker** (recommended for consistent builds):
     ```bash
     docker pull texlive/texlive
     ```

### Building the PDF

1. **Navigate to the docs directory:**
   ```bash
   cd docs/
   ```

2. **Install dependencies** (first time only):
   ```bash
   julia --project -e 'using Pkg; Pkg.instantiate()'
   ```

3. **Run the PDF build script:**
   ```bash
   julia --project make_pdf.jl
   ```

4. **Compile the LaTeX to PDF:**
   
   After running `make_pdf.jl`, LaTeX source files will be generated in `build/`. To compile to PDF:

   ```bash
   cd build/
   pdflatex JutulDarcy.jl.tex
   pdflatex JutulDarcy.jl.tex  # Run twice for references
   ```

   Or if using Docker:
   ```bash
   cd build/
   docker run --rm -v $(pwd):/work texlive/texlive pdflatex JutulDarcy.jl.tex
   ```

5. **Find your PDF:**
   The compiled PDF will be at `build/JutulDarcy.jl.pdf`

## Building Web Documentation

To build the regular web-based documentation (HTML/Vitepress):

```bash
julia --project make.jl
```

## Files

- `make.jl` - Main documentation build script (generates HTML/Vitepress documentation)
- `make_pdf.jl` - PDF documentation build script (generates LaTeX/PDF documentation using Documenter's built-in LaTeX writer)
- `Project.toml` - Julia dependencies for documentation building
- `package.json` - Node.js dependencies for Vitepress
- `src/` - Documentation source files (Markdown)

## Notes

- The PDF build includes the core documentation but **excludes dynamically generated examples** from Literate.jl to keep the PDF focused and manageable
- If you encounter LaTeX compilation errors, ensure you have a complete LaTeX distribution installed
- The `docker` platform option in `make_pdf.jl` provides the most consistent results across different systems
- For automated builds, consider using latexmk: `latexmk -pdf JutulDarcy.jl.tex`

## Troubleshooting

### Missing LaTeX packages
If you get errors about missing LaTeX packages, install them using your LaTeX distribution's package manager:
- For TeX Live: `tlmgr install <package-name>`
- For MiKTeX: Packages are usually installed automatically when needed

### Out of memory errors
If the PDF compilation runs out of memory, you may need to:
1. Reduce the documentation content in `make_pdf.jl`
2. Increase LaTeX's memory limits
3. Use a more powerful machine or cloud build service
