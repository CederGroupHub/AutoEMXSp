from PIL import Image, TiffTags

# Path to your multi-page TIFF
tiff_path = "example_particle_image.tif"

# Open the TIFF
im = Image.open(tiff_path)

# --- Extract description from page 0 ---
im.seek(0)  # Go to first image
tags = im.tag_v2
desc_str = tags.get(270, "")  # 270 = ImageDescription in TIFF spec

print("Extracted description from page 0:", repr(desc_str))

# --- Go to second image ---
im.seek(1)
second_image = im.convert("RGB")

# Save the second page using description from first page
save_path = "example2.tif"
second_image.save(
    save_path,
    format="TIFF",
    description=desc_str,
    compression=None
)

print(f"Second image saved to: {save_path}")