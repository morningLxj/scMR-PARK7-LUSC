import os
from PIL import Image, ImageEnhance

def optimize_image(image_path):
    if not os.path.exists(image_path):
        print(f"Error: File not found at {image_path}")
        return

    try:
        # Open the image
        img = Image.open(image_path)
        print(f"Original image size: {img.size}")
        print(f"Original image mode: {img.mode}")

        # Convert to RGB if necessary (for PDF saving and consistency)
        if img.mode != 'RGB':
            img = img.convert('RGB')

        # Enhance Contrast
        enhancer = ImageEnhance.Contrast(img)
        img_contrast = enhancer.enhance(1.2)  # Increase contrast by 20%

        # Enhance Sharpness
        enhancer = ImageEnhance.Sharpness(img_contrast)
        img_sharp = enhancer.enhance(1.5)  # Increase sharpness

        # Define output paths
        dir_name = os.path.dirname(image_path)
        file_name = os.path.basename(image_path)
        base_name, _ = os.path.splitext(file_name)
        
        output_png = os.path.join(dir_name, f"{base_name}_optimized.png")
        output_pdf = os.path.join(dir_name, f"{base_name}_optimized.pdf")

        # Save as 600 DPI PNG
        # DPI is metadata; for actual higher resolution, we might need to resize, 
        # but "clarity" usually implies processing, and 600dpi request usually refers to print quality metadata.
        # If the user wants physical upsizing, we would use img.resize(), but that doesn't add real detail.
        # We will stick to metadata DPI + enhancement.
        img_sharp.save(output_png, dpi=(600, 600))
        print(f"Saved optimized PNG to: {output_png}")

        # Save as PDF
        # PDF saving in PIL respects the resolution if passed or implied by image size/dpi
        img_sharp.save(output_pdf, "PDF", resolution=600.0)
        print(f"Saved optimized PDF to: {output_pdf}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # Path found from previous step
    image_path = r"d:\2026YJ\论文写作\主图\1.png"
    optimize_image(image_path)
