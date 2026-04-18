from PIL import Image, ImageDraw, ImageFont
import argparse
import sys
import os

BG_COLOR = (255, 255, 255)
LINE_COLOR = (0, 0, 0)
FILL_COLORS = [
    (200, 220, 255),
    (220, 255, 220),
    (255, 220, 220),
    (240, 240, 200),
]

MARGIN = 40
TITLE_HEIGHT = 60
TITLE_FONT_SIZE = 28
RESOLUTION_SCALE = 3

def parse_points(filename):
    if not os.path.exists(filename):
        print(f"Error: {filename} not found")
        sys.exit(1)

    with open(filename, 'r') as f:
        lines = [l.strip() for l in f.readlines() if l.strip()]

    n = int(lines[0])
    points = []

    for line in lines[1:n+1]:
        x, y = map(float, line.split())
        points.append((x, y))

    return points


def parse_triangles(filename):
    if not os.path.exists(filename):
        print(f"Error: {filename} not found")
        sys.exit(1)

    triangles = []

    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            i, j, k = map(int, parts)
            triangles.append((i, j, k))

    return triangles


def compute_bounds(points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return min(xs), max(xs), min(ys), max(ys)


def draw_mesh(points, triangles, output_file, scale=5):
    min_x, max_x, min_y, max_y = compute_bounds(points)

    width = int((max_x - min_x) * scale * RESOLUTION_SCALE) + 2 * MARGIN * RESOLUTION_SCALE
    height = int((max_y - min_y) * scale * RESOLUTION_SCALE) + 2 * MARGIN * RESOLUTION_SCALE + TITLE_HEIGHT * RESOLUTION_SCALE

    img = Image.new("RGB", (width, height), BG_COLOR)
    draw = ImageDraw.Draw(img)

    title = f"{len(points)} points, {len(triangles)} triangles"

    try:
        font = ImageFont.truetype("DejaVuSans.ttf", TITLE_FONT_SIZE * RESOLUTION_SCALE)
    except:
        font = ImageFont.load_default()

    draw.text((MARGIN, 10), title, fill=(0, 0, 0), font=font)

    def transform(p):
        x, y = p
        tx = (x - min_x) * scale * RESOLUTION_SCALE + MARGIN * RESOLUTION_SCALE
        ty = (y - min_y) * scale * RESOLUTION_SCALE + MARGIN * RESOLUTION_SCALE + TITLE_HEIGHT * RESOLUTION_SCALE
        return (tx, ty)

    #* draw triangles
    for idx, (i, j, k) in enumerate(triangles):
        poly = [
            transform(points[i]),
            transform(points[j]),
            transform(points[k])
        ]

        color = FILL_COLORS[idx % len(FILL_COLORS)]

        draw.polygon(poly, fill=color)
        draw.line([poly[0], poly[1], poly[2], poly[0]], fill=LINE_COLOR, width=1*RESOLUTION_SCALE)

    #* draw points
    r = 2 * RESOLUTION_SCALE
    for p in points:
        x, y = transform(p)
        draw.ellipse((x-r, y-r, x+r, y+r), fill=(0, 0, 0))

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    img.save(output_file)
    print(f"Saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--points', type=str, default='inputs/dummy.txt')
    parser.add_argument('--triangles', type=str, default='outputs/dummy.txt')
    parser.add_argument('--out', type=str, default='outputs/visuals/mesh.png')
    parser.add_argument('--scale', type=int, default=6)
    args = parser.parse_args()

    points = parse_points(args.points)
    triangles = parse_triangles(args.triangles)

    draw_mesh(points, triangles, args.out, args.scale)