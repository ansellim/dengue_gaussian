#!/usr/bin/env bash
# Split the combined slides+notes PDF (produced by `\setbeameroption{show notes}`)
# into a slides-only PDF and a speaker-notes-only PDF.
#
# Scans slides.tex to figure out which frames have a `\note{...}` following them.
# Frames without a note contribute only a slide page; frames with a note contribute
# a slide page immediately followed by a note page. This handles the case where
# some slides (e.g. the Key References slide) have no speaker notes.
set -euo pipefail

cd "$(dirname "$0")"

INPUT="${1:-slides.pdf}"
TEX="${INPUT%.pdf}.tex"
SLIDES_OUT="slides_only.pdf"
NOTES_OUT="speaker_notes_only.pdf"

if [[ ! -f "$INPUT" ]]; then
  echo "error: $INPUT not found. Render slides.qmd first." >&2
  exit 1
fi
if [[ ! -f "$TEX" ]]; then
  echo "error: $TEX not found. Quarto beamer output must keep the .tex (keep-tex: true)." >&2
  exit 1
fi

# Walk slides.tex in order. For every frame emit its slide page; if a \note{ block
# appears before the next frame, emit the following page as a note page.
mapfile -t pages < <(
  awk '
    BEGIN { page = 0; pending_note = 0 }
    function flush_slide() {
      if (pending_frame) {
        page += 1
        print "S " page
        pending_note = 1
        pending_frame = 0
      }
    }
    /\\frame\{\\titlepage\}/ {
      flush_slide()
      pending_frame = 1
      pending_note_allowed = 0  # title has no \note
      next
    }
    /\\begin\{frame\}/ {
      flush_slide()
      pending_frame = 1
      pending_note_allowed = 1
      next
    }
    /^\\note\{/ {
      if (pending_frame) {
        page += 1
        print "S " page
        pending_frame = 0
      }
      if (pending_note_allowed) {
        page += 1
        print "N " page
        pending_note_allowed = 0
      }
      next
    }
    END {
      flush_slide()
    }
  ' "$TEX"
)

slide_pages=()
note_pages=()
for entry in "${pages[@]}"; do
  kind="${entry%% *}"
  num="${entry##* }"
  if [[ $kind == S ]]; then
    slide_pages+=("$num")
  else
    note_pages+=("$num")
  fi
done

total_pdf=$(pdfinfo "$INPUT" | awk '/^Pages:/ {print $2}')
total_parsed=${#pages[@]}
if [[ "$total_parsed" -ne "$total_pdf" ]]; then
  echo "warning: parsed $total_parsed pages from $TEX but $INPUT has $total_pdf pages." >&2
  echo "         the .tex and .pdf may be out of sync; re-render slides.qmd." >&2
fi

if [[ ${#slide_pages[@]} -eq 0 ]]; then
  echo "error: found no slide pages in $TEX." >&2
  exit 1
fi

pdftk "$INPUT" cat "${slide_pages[@]}" output "$SLIDES_OUT"
if [[ ${#note_pages[@]} -gt 0 ]]; then
  pdftk "$INPUT" cat "${note_pages[@]}" output "$NOTES_OUT"
else
  echo "warning: no note pages found; $NOTES_OUT not written." >&2
  rm -f "$NOTES_OUT"
fi

slides_n=$(pdfinfo "$SLIDES_OUT" | awk '/^Pages:/ {print $2}')
echo "Wrote:"
echo "  $SLIDES_OUT ($slides_n pages)"
if [[ -f "$NOTES_OUT" ]]; then
  notes_n=$(pdfinfo "$NOTES_OUT" | awk '/^Pages:/ {print $2}')
  echo "  $NOTES_OUT  ($notes_n pages)"
fi
