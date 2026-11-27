from PIL import Image

def merge_2x2(image_paths, output_path, padding=20, bg_color="white"):
    imgs = [Image.open(p) for p in image_paths]

    # 统一尺寸：按四张图中的最大宽高补白
    max_w = max(img.width for img in imgs)
    max_h = max(img.height for img in imgs)

    norm_imgs = []
    for img in imgs:
        canvas = Image.new("RGB", (max_w, max_h), bg_color)
        x = (max_w - img.width) // 2
        y = (max_h - img.height) // 2
        canvas.paste(img, (x, y))
        norm_imgs.append(canvas)

    panel_w = max_w * 2 + padding * 3
    panel_h = max_h * 2 + padding * 3
    panel = Image.new("RGB", (panel_w, panel_h), bg_color)

    positions = [
        (padding, padding),
        (padding * 2 + max_w, padding),
        (padding, padding * 2 + max_h),
        (padding * 2 + max_w, padding * 2 + max_h),
    ]

    for img, pos in zip(norm_imgs, positions):
        panel.paste(img, pos)

    panel.save(output_path)

merge_2x2(
    [
        "106788.deeprm.beta_hek_vs_beta_hela.png",
        "106788.dorado.beta_hek_vs_beta_hela.png",
        "106788.m6anet.beta_hek_vs_beta_hela.png",
        "106788.singlemod.beta_hek_vs_beta_hela.png",
    ],
    "106788.beta_hek_vs_beta_hela.panel.png"
)

merge_2x2(
    [
        "80229.deeprm.beta_hek_vs_beta_hela.png",
        "80229.dorado.beta_hek_vs_beta_hela.png",
        "80229.m6anet.beta_hek_vs_beta_hela.png",
        "80229.singlemod.beta_hek_vs_beta_hela.png",
    ],
    "80229.beta_hek_vs_beta_hela.panel.png"
)
