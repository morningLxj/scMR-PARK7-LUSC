import os
import fitz

base = "My_MR_Project/TCGA_Results"
files = [
    os.path.join(base, "TCGA-LUAD_TMEM50A_Survival.pdf"),
    os.path.join(base, "TCGA-LUSC_CTSW_Survival.pdf"),
    os.path.join(base, "TCGA-LUAD_PARK7_Survival.pdf"),
]
docs = [fitz.open(p) for p in files]
# use page 1 (second page) when available, else page 0
page_indices = [1 if d.page_count >= 2 else 0 for d in docs]
pages = [d[idx] for d, idx in zip(docs, page_indices)]
sizes = [(p.rect.width, p.rect.height) for p in pages]
target_h = max(h for (_, h) in sizes)
scales = [target_h / h for (_, h) in sizes]
scaled_ws = [w * s for (w, _), s in zip(sizes, scales)]
gap = 18.0
margin = 18.0
total_w = sum(scaled_ws) + margin * 2 + gap * (len(pages) - 1)
total_h = target_h + margin * 2
out_pdf = "My_MR_Project/Figure5_TCGA_Survival_Composite.pdf"
out_png = "My_MR_Project/Figure5_TCGA_Survival_Composite_600dpi.png"

doc_out = fitz.open()
page_out = doc_out.new_page(width=total_w, height=total_h)
x = margin
labels = ["A", "B", "C"]
for i, (d, w, s, idx) in enumerate(zip(docs, scaled_ws, scales, page_indices)):
    rect = fitz.Rect(x, margin, x + w, margin + target_h)
    page_out.show_pdf_page(rect, d, idx)
    page_out.insert_text((x + 8, margin - 6), labels[i], fontsize=16, fontname="helv")
    x += w + gap
doc_out.save(out_pdf)
doc_out.close()
for d in docs:
    d.close()

doc = fitz.open(out_pdf)
pg = doc[0]
zoom = 600.0 / 72.0
mat = fitz.Matrix(zoom, zoom)
pix = pg.get_pixmap(matrix=mat, alpha=False)
pix.save(out_png)
doc.close()
print(f"Composite saved: {out_pdf} and {out_png}")

