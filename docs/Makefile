.PHONY: help pdf clean

help:
	@echo "JutulDarcy.jl Documentation Build"
	@echo "=================================="
	@echo ""
	@echo "Available targets:"
	@echo "  make pdf    - Build PDF documentation using Pandoc"
	@echo "  make clean  - Remove generated PDF documentation"
	@echo "  make help   - Show this help message"
	@echo ""
	@echo "For more information, see README_PDF.md"

pdf:
	@echo "Building PDF documentation..."
	./build_pdf_pandoc.sh

clean:
	@echo "Cleaning up generated PDF..."
	rm -f JutulDarcy_Documentation.pdf
	@echo "Done."
