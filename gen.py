import argparse
import math
import random


def generate_points(n, evenness):
    """
    evenness:
        0.0 -> perfect grid
        1.0 -> fully random
    """

    points = []

    if evenness >= 1.0:
        for _ in range(n):
            x = random.random()
            y = random.random()
            points.append((x, y))
        return points

    # grid-based generation
    grid_size = math.ceil(math.sqrt(n))
    cell_size = 1000.0 / grid_size

    count = 0
    for i in range(grid_size):
        for j in range(grid_size):
            if count >= n:
                break

            base_x = (i + 0.5) * cell_size
            base_y = (j + 0.5) * cell_size

            jitter = evenness * cell_size * 0.5

            x = base_x + random.uniform(-jitter, jitter)
            y = base_y + random.uniform(-jitter, jitter)

            # clamp
            x = min(max(x, 0.0), 1.0)
            y = min(max(y, 0.0), 1.0)

            points.append((x, y))
            count += 1

        if count >= n:
            break

    return points


def write_file(filename, points):
    with open(filename, "w") as f:
        f.write(f"{len(points)}\n")
        for x, y in points:
            f.write(f"{x:.6f} {y:.6f}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--evenness", type=float, default=0.5)
    parser.add_argument("--out", type=str, default="points.txt")

    args = parser.parse_args()

    points = generate_points(args.n, args.evenness)
    write_file(args.out, points)

    print(f"Wrote {len(points)} points to {args.out}")


if __name__ == "__main__":
    main()