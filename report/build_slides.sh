#!/usr/bin/env bash
# Render slides.qmd to a combined slides+notes PDF, then split it into
# slides-only and speaker-notes-only PDFs.
set -euo pipefail

cd "$(dirname "$0")"

QUARTO="${QUARTO:-quarto}"
if ! command -v "$QUARTO" >/dev/null 2>&1; then
  QUARTO="/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto"
fi

echo "==> Rendering slides.qmd"
"$QUARTO" render slides.qmd

echo "==> Splitting slides.pdf into slides-only and notes-only"
./split_slides_notes.sh slides.pdf

echo "==> Done"